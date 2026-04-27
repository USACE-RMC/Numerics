using System;
using System.Collections.Generic;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Inverse Gamma distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Inverse-gamma_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class InverseGamma : UnivariateDistributionBase
    {
        /// <summary>
        /// Constructs an Inverse-Gamma distribution with scale β = 0.5 and shape α = 2.
        /// </summary>
        public InverseGamma()
        {
            SetParameters([0.5d, 2d]);
        }

        /// <summary>
        /// Constructs a new Inverse-Gamma distribution with a given scale and shape.
        /// </summary>
        /// <param name="scale">The scale parameter β (beta).</param>
        /// <param name="shape">The shape parameter α (alpha).</param>
        public InverseGamma(double scale, double shape)
        {
            SetParameters([scale, shape]);
        }

        private double _beta;
        private double _alpha;

        /// <summary>
        /// Gets and sets the scale parameter β (beta).
        /// </summary>
        public double Beta
        {
            get { return _beta; }
            set
            {
                _parametersValid = ValidateParameters([value, Alpha], false) is null;
                _beta = value;
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
                _parametersValid = ValidateParameters([Beta, value], false) is null;
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
            get { return UnivariateDistributionType.InverseGamma; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Inverse-Gamma"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "Inv-G"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Scale (β)";
                parmString[1, 0] = "Shape (α)";
                parmString[0, 1] = Beta.ToString();
                parmString[1, 1] = Alpha.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["β", "α"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Beta), nameof(Alpha)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Beta, Alpha]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get 
            {
                if (Alpha <= 1)
                    return double.NaN;
                return Beta / (Alpha - 1d); 
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
            get { return Beta / (Alpha + 1d); }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                if (Alpha <= 2.0d)
                    return double.NaN;
                return Beta / (Math.Abs(Alpha - 1.0d) * Math.Sqrt(Alpha - 2.0d));
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {
                if (Alpha <= 3.0d)
                    return double.NaN;
                return 4.0d * Math.Sqrt(Alpha - 2.0d) / (Alpha - 3.0d);
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                if (Alpha <= 4.0d)
                    return double.NaN;
                return 3.0d + 6.0d * (5.0d * Alpha - 11.0d) / ((Alpha - 3.0d) * (Alpha - 4.0d));
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
            get { return [0.0d, 0.0d]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity, double.PositiveInfinity]; }
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            Beta = parameters[0];
            Alpha = parameters[1];
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (double.IsNaN(parameters[0]) || double.IsInfinity(parameters[0]) || parameters[0] <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Beta), "The scale parameter β (beta)) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Beta), "The scale parameter β (beta) must be positive.");
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
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([Beta, Alpha], true);
            if (x < Minimum || x > Maximum) return 0.0d;
            return Math.Pow(Beta, Alpha) / Gamma.Function(Alpha) * Math.Pow(x, -Alpha - 1d) * Math.Exp(-Beta / x);
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([Beta, Alpha], true);
            if (x <= Minimum)
                return 0d;
            if (x >= Maximum)
                return 1d;
            return Gamma.UpperIncomplete(Alpha, Beta / x);
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (probability < 0.0d || probability > 1.0d)
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            if (probability == 0.0d)
                return Minimum;
            if (probability == 1.0d)
                return Maximum;
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([Beta, Alpha], true);
            return Beta / Gamma.InverseUpperIncomplete(Alpha, probability);
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new InverseGamma(Beta, Alpha);
        }

    }
}