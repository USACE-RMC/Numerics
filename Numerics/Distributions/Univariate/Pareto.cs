using System;
using System.Collections.Generic;
using Numerics.Data.Statistics;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Pareto distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Pareto_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class Pareto : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IBootstrappable
    {
  
        /// <summary>
        /// Constructs a Pareto distribution with scale = 1 and shape = 10.
        /// </summary>
        public Pareto()
        {
            SetParameters([1d, 10d]);
        }

        /// <summary>
        /// Constructs a Pareto distribution with the given parameters Xm and α.
        /// </summary>
        /// <param name="scale">The scale parameter Xm.</param>
        /// <param name="shape">The shape parameter α (alpha).</param>
        public Pareto(double scale, double shape)
        {
            SetParameters([scale, shape]);
        }

        private double _Xm;
        private double _alpha;

        /// <summary>
        /// Gets and sets the scale parameter Xm.
        /// </summary>
        public double Xm
        {
            get { return _Xm; }
            set
            {
                _parametersValid = ValidateParameters(new[] { value, _alpha }, false) is null;
                _Xm = value;
            }
        }

        /// <summary>
        /// Gets and sets the shape parameter α (alpha).
        /// </summary>
        public double Alpha
        {
            get { return _alpha; }
            set
            {
                _parametersValid = ValidateParameters(new[] { Xm, value }, false) is null;
                _alpha = value;
            }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 2; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.Pareto; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Pareto"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "PA"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Scale (Xm)";
                parmString[1, 0] = "Shape (α)";
                parmString[0, 1] = Xm.ToString();
                parmString[1, 1] = Alpha.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["Xm", "α"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Xm), nameof(Alpha)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Xm, Alpha]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                if (Alpha <= 1d)
                    return double.PositiveInfinity;
                return Alpha * Xm / (Alpha - 1.0d);
            }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return Xm * Math.Pow(2.0d, 1.0d / Alpha); }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get { return Xm; }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                if (Alpha <= 2d)
                    return double.PositiveInfinity;
                return Xm * Math.Sqrt(Alpha) / (Math.Abs(Alpha - 1.0d) * Math.Sqrt(Alpha - 2.0d));
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {
                if (Alpha <= 3d)
                    return double.NaN;
                return 2.0d * (Alpha + 1.0d) / (Alpha - 3.0d) * Math.Sqrt((Alpha - 2.0d) / Alpha);
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                if (Alpha <= 4d)
                    return double.NaN;
                double num = 6d * (Alpha * Alpha * Alpha + Alpha * Alpha - 6d * Alpha - 2d);
                double den = Alpha * (Alpha - 3d) * (Alpha - 4d);
                return 3.0d + num / den;
            }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return Xm; }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return double.PositiveInfinity; }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [0.0d, 0.0d]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity, double.PositiveInfinity]; }
        }

        /// <inheritdoc/>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample has fewer than two observations, contains a
        /// nonfinite or nonpositive value, or is constant.</exception>
        /// <exception cref="NotImplementedException"><paramref name="estimationMethod"/> is not
        /// <see cref="ParameterEstimationMethod.MaximumLikelihood"/>.</exception>
        public void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateSample(sample, 2, true);
            if (estimationMethod == ParameterEstimationMethod.MaximumLikelihood)
            {
                SetParameters(MLE(sample));
            }
            else
            {
                throw new NotImplementedException();
            }
        }

        /// <inheritdoc/>
        public IUnivariateDistribution Bootstrap(ParameterEstimationMethod estimationMethod, int sampleSize, int seed = -1)
        {
            var newDistribution = new Pareto(Xm, Alpha);
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The initial values are the closed-form maximum likelihood estimates. The scale is bounded above by the
        /// sample minimum, because every observation must lie in the support [Xm, ∞).
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample has fewer than two observations, contains a
        /// nonfinite or nonpositive value, or is constant.</exception>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            var initialVals = MLE(sample);
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];
            DistributionNumerics.PositiveParameterBounds(initialVals[0], out lowerVals[0], out _);
            upperVals[0] = initialVals[0];
            DistributionNumerics.PositiveParameterBounds(initialVals[1], out lowerVals[1], out upperVals[1]);
            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The maximum likelihood estimates have the closed form Xm = min(x) and α = n / Σ ln(xᵢ/Xm)
        /// (Johnson, Kotz and Balakrishnan, 1994, Continuous Univariate Distributions, Vol. 1, Ch. 20); they equal
        /// scipy.stats.pareto.fit with the location fixed at zero.
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample has fewer than two observations, contains a
        /// nonfinite or nonpositive value, or is constant.</exception>
        public double[] MLE(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 2, true);
            double scale = Statistics.Minimum(sample);
            double sum = 0d;
            for (int i = 0; i < sample.Count; i++)
                sum += Math.Log(sample[i] / scale);
            return [scale, sample.Count / sum];
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            Xm = parameters[0];
            Alpha = parameters[1];
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (double.IsNaN(parameters[0]) || double.IsInfinity(parameters[0]) || parameters[0] <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Xm), "The scale parameter Xm must be positive.");
                return new ArgumentOutOfRangeException(nameof(Xm), "The scale parameter Xm must be positive.");
            }

            if (double.IsNaN(parameters[1]) || double.IsInfinity(parameters[1]) || parameters[1] <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Alpha), "The shape parameter α (alpha) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Alpha), "The shape parameter α (alpha) must be positive.");
            }
            return null!;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated as α/x · (Xm/x)^α, so the density stays finite when Xm^α or x^(α+1) alone would overflow.
        /// </remarks>
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(new[] { Xm, Alpha }, true);
            if (x < Minimum || x > Maximum) return 0.0d;
            return Alpha / x * Math.Pow(Xm / x, Alpha);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated as ln(α) − α·ln(x/Xm) − ln(x), which stays finite after the density underflows.
        /// </remarks>
        public override double LogPDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(new[] { Xm, Alpha }, true);
            if (x < Minimum || x > Maximum) return double.NegativeInfinity;
            return Math.Log(Alpha) - Alpha * LogScaledExcess(x) - Math.Log(x);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated as −expm1(−α·ln(x/Xm)) with ln(x/Xm) = log1p((x − Xm)/Xm) near Xm, so the
        /// probability keeps full relative precision just above the scale parameter.
        /// </remarks>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(new[] { Xm, Alpha }, true);
            if (x <= Minimum)
                return 0d;
            if (x >= Maximum)
                return 1d;
            return -Tools.Expm1(-Alpha * LogScaledExcess(x));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated directly as (Xm/x)^α = exp(−α·ln(x/Xm)), so far upper-tail probabilities are not lost
        /// to the subtraction 1 − CDF.
        /// </remarks>
        public override double CCDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(new[] { Xm, Alpha }, true);
            if (x <= Minimum)
                return 1d;
            return Math.Exp(-Alpha * LogScaledExcess(x));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated as −α·ln(x/Xm), which stays finite after the survival probability underflows.
        /// </remarks>
        public override double LogCCDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(new[] { Xm, Alpha }, true);
            if (x <= Minimum)
                return 0d;
            return -Alpha * LogScaledExcess(x);
        }

        /// <summary>
        /// Computes ln(x/Xm) for x at or above the scale parameter.
        /// </summary>
        /// <param name="x">A value at or above <see cref="Xm"/>.</param>
        /// <returns>The natural logarithm of x/Xm.</returns>
        /// <remarks>
        /// Near Xm the logarithm is formed as log1p((x − Xm)/Xm), where the difference is exact, so it keeps
        /// full relative precision; farther out the ratio is used directly, and when the ratio overflows the
        /// difference of logarithms is used.
        /// </remarks>
        private double LogScaledExcess(double x)
        {
            double excess = (x - Xm) / Xm;
            if (excess <= 1d)
                return Tools.Log1p(excess);
            double ratio = x / Xm;
            return double.IsPositiveInfinity(ratio) ? Math.Log(x) - Math.Log(Xm) : Math.Log(ratio);
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (probability < 0.0d || probability > 1.0d)
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(new[] { Xm, Alpha }, true);
            if (probability == 0.0d)
                return Minimum;
            if (probability == 1.0d)
                return Maximum;
            return Xm * Math.Pow(1.0d - probability, -1.0d / Alpha);
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new Pareto(Xm, Alpha);
        }
          
    }
}