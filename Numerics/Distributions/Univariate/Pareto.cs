﻿using System;
using System.Collections.Generic;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Pareto distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     Authors:
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Pareto_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class Pareto : UnivariateDistributionBase
    {
  
        /// <summary>
        /// Constructs a Pareto distribution with scale = 1 and shape = 10.
        /// </summary>
        public Pareto()
        {
            SetParameters(new[] { 1d, 10d });
        }

        /// <summary>
        /// Constructs a Pareto distribution with the given parameters Xm and α.
        /// </summary>
        /// <param name="scale">The scale parameter Xm.</param>
        /// <param name="shape">The shape parameter α (alpha).</param>
        public Pareto(double scale, double shape)
        {
            SetParameters(new[] { scale, shape });
        }

        private bool _parametersValid = true;
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

        /// <summary>
        /// Returns the number of distribution parameters.
        /// </summary>
        public override int NumberOfParameters
        {
            get { return 2; }
        }

        /// <summary>
        /// Returns the univariate distribution type.
        /// </summary>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.Pareto; }
        }

        /// <summary>
        /// Returns the name of the distribution type as a string.
        /// </summary>
        public override string DisplayName
        {
            get { return "Pareto"; }
        }

        /// <summary>
        /// Returns the short display name of the distribution as a string.
        /// </summary>
        public override string ShortDisplayName
        {
            get { return "PA"; }
        }

        /// <summary>
        /// Get distribution parameters in 2-column array of string.
        /// </summary>
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

        /// <summary>
        /// Gets the short form parameter names.
        /// </summary>
        public override string[] ParameterNamesShortForm
        {
            get { return new[] { "Xm", "α" }; }
        }

        /// <summary>
        /// Gets the full parameter names.
        /// </summary>
        public override string[] GetParameterPropertyNames
        {
            get { return new[] { nameof(Xm), nameof(Alpha) }; }
        }

        /// <summary>
        /// Get an array of parameters.
        /// </summary>
        public override double[] GetParameters
        {
            get { return new[] { Xm, Alpha }; }
        }

        /// <summary>
        /// Determines whether the parameters are valid or not.
        /// </summary>
        public override bool ParametersValid
        {
            get { return _parametersValid; }
        }

        /// <summary>
        /// Gets the mean of the distribution.
        /// </summary>
        public override double Mean
        {
            get
            {
                if (Alpha <= 1d)
                    return double.PositiveInfinity;
                return Alpha * Xm / (Alpha - 1.0d);
            }
        }

        /// <summary>
        /// Gets the median of the distribution.
        /// </summary>
        public override double Median
        {
            get { return Xm * Math.Pow(2.0d, 1.0d / Alpha); }
        }

        /// <summary>
        /// Gets the mode of the distribution.
        /// </summary>
        public override double Mode
        {
            get { return Xm; }
        }

        /// <summary>
        /// Gets the standard deviation of the distribution.
        /// </summary>
        public override double StandardDeviation
        {
            get
            {
                if (Alpha <= 2d)
                    return double.PositiveInfinity;
                return Xm * Math.Sqrt(Alpha) / (Math.Abs(Alpha - 1.0d) * Math.Sqrt(Alpha - 2.0d));
            }
        }

        /// <summary>
        /// Gets the skew of the distribution.
        /// </summary>
        public override double Skew
        {
            get
            {
                if (Alpha <= 3d)
                    return double.NaN;
                return 2.0d * (Alpha + 1.0d) / (Alpha - 3.0d) * Math.Sqrt((Alpha - 2.0d) / Alpha);
            }
        }

        /// <summary>
        /// Gets the kurtosis of the distribution.
        /// </summary>
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

        /// <summary>
        /// Gets the minimum of the distribution.
        /// </summary>
        public override double Minimum
        {
            get { return Xm; }
        }

        /// <summary>
        /// Gets the maximum of the distribution.
        /// </summary>
        public override double Maximum
        {
            get { return double.PositiveInfinity; }
        }

        /// <summary>
        /// Gets the minimum values allowable for each parameter.
        /// </summary>
        public override double[] MinimumOfParameters
        {
            get { return new[] { 0.0d, 0.0d }; }
        }

        /// <summary>
        /// Gets the maximum values allowable for each parameter.
        /// </summary>
        public override double[] MaximumOfParameters
        {
            get { return new[] { double.PositiveInfinity, double.PositiveInfinity }; }
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="parameters">A list of parameters.</param>
        public override void SetParameters(IList<double> parameters)
        {
            Xm = parameters[0];
            Alpha = parameters[1];
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="parameters">A list of parameters.</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (parameters[0] <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Xm), "The scale parameter Xm must be positive.");
                return new ArgumentOutOfRangeException(nameof(Xm), "The scale parameter Xm must be positive.");
            }

            if (parameters[1] <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Alpha), "The shape parameter α (alpha) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Alpha), "The shape parameter α (alpha) must be positive.");
            }
            return null;
        }
     
        /// <summary>
        /// Gets the Probability Density Function (PDF) of the distribution evaluated at a point x.
        /// </summary>
        /// <param name="x">A single point in the distribution range.</param>
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(new[] { Xm, Alpha }, true);
            if (x < Minimum || x > Maximum) return 0.0d;
            return Alpha * Math.Pow(Xm, Alpha) / Math.Pow(x, Alpha + 1d);
        }

        /// <summary>
        /// Gets the Cumulative Distribution Function (CDF) for the distribution evaluated at a point x.
        /// </summary>
        /// <param name="x">A single point in the distribution range.</param>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(new[] { Xm, Alpha }, true);
            if (x <= Minimum)
                return 0d;
            if (x >= Maximum)
                return 1d;
            return 1d - Math.Pow(Xm / x, Alpha);
        }

        /// <summary>
        /// Gets the Inverse Cumulative Distribution Function (ICFD) of the distribution evaluated at a probability.
        /// </summary>
        /// <param name="probability">Probability between 0 and 1.</param>
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
                ValidateParameters(new[] { Xm, Alpha }, true);
            return Xm * Math.Pow(1.0d - probability, -1.0d / Alpha);
        }
     
        /// <summary>
        /// Creates a copy of the distribution.
        /// </summary>
        public override UnivariateDistributionBase Clone()
        {
            return new Pareto(Xm, Alpha);
        }
          
    }
}