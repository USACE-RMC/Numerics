using System;
using System.Collections.Generic;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Inverse Chi-Squared (Inv-χ²) probability distribution.
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
    /// <see href = "https://en.wikipedia.org/wiki/Inverse-chi-squared_distribution" />
    /// <see href = "https://en.wikipedia.org/wiki/Scaled_inverse_chi-squared_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class InverseChiSquared : UnivariateDistributionBase
    {
   
        /// <summary>
        /// Constructs an Inverse Chi-Squared distribution with 10 degrees of freedom.
        /// </summary>
        public InverseChiSquared()
        {
            SetParameters([10d, 1d]);
        }

        /// <summary>
        /// Constructs an Inverse Chi-Squared distribution with given degrees of freedom.
        /// </summary>
        /// <param name="degreesOfFreedom">The degrees of freedom for the distribution.</param>
        /// <param name="scale">The scale parameter τ² of the scaled inverse chi-squared distribution (the <see cref="Sigma"/> property).</param>
        public InverseChiSquared(int degreesOfFreedom, double scale)
        {
            SetParameters([degreesOfFreedom, scale]);
        }

        private int _degreesOfFreedom;
        private double _sigma;

        /// <summary>
        /// Gets and sets the degrees of freedom ν (nu) of the distribution.
        /// </summary>
        public int DegreesOfFreedom
        {
            get { return _degreesOfFreedom; }
            set
            {
                _parametersValid = ValidateParameters([value, Sigma], false) is null;
                _degreesOfFreedom = value;
            }
        }

        /// <summary>
        /// Gets and sets the scale parameter τ², stored in the property named Sigma.
        /// </summary>
        /// <remarks>
        /// This value is the scale parameter τ² of the scaled inverse chi-squared distribution, which is the
        /// inverse-gamma distribution with shape ν/2 and scale ν·τ²/2. The mean is ν·τ²/(ν - 2) for ν &gt; 2.
        /// </remarks>
        public double Sigma
        {
            get { return _sigma; }
            set
            {
                _parametersValid = ValidateParameters([DegreesOfFreedom, value], false) is null;
                _sigma = value;
            }
        }

        /// <summary>
        /// Gets the value 1/σ², where σ is <see cref="Sigma"/>.
        /// </summary>
        /// <remarks>
        /// Because <see cref="Sigma"/> is the scale τ², this value equals 1/τ⁴, not the reciprocal of the scale.
        /// </remarks>
        public double Precision
        {
            get { return 1d / (Sigma * Sigma); }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 2; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.InverseChiSquared; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Inverse Chi-Squared (χ²)"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "Inv-χ²"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Degrees of Freedom (ν)";
                parmString[1, 0] = "Scale (σ)";
                parmString[0, 1] = DegreesOfFreedom.ToString();
                parmString[1, 1] = Sigma.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["ν", "σ"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames => [nameof(DegreesOfFreedom), nameof(Sigma)];

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [DegreesOfFreedom, Sigma]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                if (DegreesOfFreedom <= 2)
                    return double.NaN;
                return DegreesOfFreedom  * Sigma / (DegreesOfFreedom - 2);
            }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return InverseCDF(0.5d); }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get { return DegreesOfFreedom * Sigma / (DegreesOfFreedom + 2); }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                if (DegreesOfFreedom <= 4)
                    return double.NaN;
                int v = DegreesOfFreedom;
                double t2 = Sigma;
                return Math.Sqrt(2.0d * v * v * t2 * t2 / ((v - 2) * (v - 2) * (v - 4)));
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {
                if (DegreesOfFreedom <= 6)
                    return double.NaN;
                int v = DegreesOfFreedom;
                return 4.0d / (v - 6) * Math.Sqrt(2.0d * (v - 4));
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                if (DegreesOfFreedom <= 8)
                    return double.NaN;
                int v = DegreesOfFreedom;
                return 3.0d + 12.0d * (5.0d * v - 22.0d) / ((v - 6) * (v - 8));
            }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return 0.0d; }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return double.PositiveInfinity; }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [1.0d, 0.0d]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity, double.PositiveInfinity]; }
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            _degreesOfFreedom = (int)parameters[0];
            _sigma = parameters[1];
            _parametersValid = ValidateParameters(parameters, false) is null;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (double.IsNaN(parameters[0]) || double.IsInfinity(parameters[0]) || parameters[0] < 1.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(DegreesOfFreedom), "The degrees of freedom ν (nu) must greater than or equal to one.");
                return new ArgumentOutOfRangeException(nameof(DegreesOfFreedom), "The degrees of freedom ν (nu) must greater than or equal to one.");
            }
            if (double.IsNaN(parameters[1]) || double.IsInfinity(parameters[1]) || parameters[1] <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Sigma), "The scale parameter σ (sigma) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Sigma), "The scale parameter σ (sigma) must be positive.");
            }
            return null!;
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            return Math.Exp(LogPDF(x));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated in log space through the unit-scale gamma density of y = ν·σ / (2x), using
        /// log f(x) = log g(y; ν/2) + log y - log x, so large degrees of freedom keep a finite density
        /// and the density is zero at x = 0.
        /// </remarks>
        public override double LogPDF(double x)
        {
            // validate parameters
            if (_parametersValid == false)
                ValidateParameters([DegreesOfFreedom, Sigma], true);
            if (x <= Minimum || double.IsPositiveInfinity(x)) return double.NegativeInfinity;
            double v = DegreesOfFreedom;
            double t2 = Sigma;
            double scale = v * t2 / 2d;
            double y = scale / x;
            if (double.IsPositiveInfinity(y)) return double.NegativeInfinity;
            // When the scaled variate underflows, exp(-y) is one and the power-law tail remains.
            if (y == 0d) return v / 2d * Math.Log(scale) - Gamma.LogGamma(v / 2d) - (v / 2d + 1d) * Math.Log(x);
            return DistributionNumerics.GammaLogDensity(v / 2d, y) + Math.Log(y) - Math.Log(x);
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([DegreesOfFreedom, Sigma], true);
            if (x <= Minimum)
                return 0d;
            double v = DegreesOfFreedom;
            double t2 = Sigma;
            return Gamma.UpperIncomplete(v / 2.0d, v * t2 / (2.0d * x));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated as the log regularized upper incomplete gamma function Q(ν/2, ν·σ / (2x)).
        /// </remarks>
        public override double LogCDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([DegreesOfFreedom, Sigma], true);
            if (x <= Minimum)
                return double.NegativeInfinity;
            if (x >= Maximum)
                return 0d;
            double v = DegreesOfFreedom;
            double y = v * Sigma / (2.0d * x);
            return y == 0d ? Tools.Log1mExp(LogCCDF(x)) : DistributionNumerics.GammaLogSurvival(v / 2.0d, y);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated as the regularized lower incomplete gamma function P(ν/2, ν·σ / (2x)), which resolves
        /// upper-tail probabilities below machine epsilon.
        /// </remarks>
        public override double CCDF(double x)
        {
            return Math.Exp(LogCCDF(x));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated as the log regularized lower incomplete gamma function P(ν/2, ν·σ / (2x)).
        /// </remarks>
        public override double LogCCDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([DegreesOfFreedom, Sigma], true);
            if (x <= Minimum)
                return 0d;
            if (x >= Maximum)
                return double.NegativeInfinity;
            double v = DegreesOfFreedom;
            double scale = v * Sigma / 2.0d;
            double y = v * Sigma / (2.0d * x);
            // When the scaled variate underflows, the leading series term P(a, y) = y^a / Γ(a + 1) is exact to double precision.
            return y == 0d ? v / 2.0d * (Math.Log(scale) - Math.Log(x)) - Gamma.LogGammaOnePlus(v / 2.0d)
                : DistributionNumerics.GammaLogCDF(v / 2.0d, y);
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (!(probability >= 0.0d && probability <= 1.0d))
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            if (probability == 0.0d)
                return Minimum;
            if (probability == 1.0d)
                return Maximum;
            // validate parameters
            if (_parametersValid == false)
                ValidateParameters([DegreesOfFreedom, Sigma], true);
            double v = DegreesOfFreedom;
            double t2 = Sigma;
            return v * t2 / (2.0d * Gamma.InverseUpperIncomplete(v / 2.0d, probability));
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new InverseChiSquared(DegreesOfFreedom, Sigma);
        }

    }
}