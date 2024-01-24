﻿using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using Numerics.Data;
using Numerics.Data.Statistics;
using Numerics.Mathematics;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Univariate Empirical distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     Authors:
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// This distribution specifies a cumulative distribution with n points. The range of the distribution
    /// is set by the minimum and maximum arguments. Each point on the cumulative curve has an X value and
    /// a probability. Points on the cumulative curve must be entered with increasing value and increasing
    /// probability. Even though the (X,p) points define the distribution, and value between the minimum
    /// and maximum can be returned.
    /// </para>
    /// <para>
    /// References:
    /// <list type="bullet">
    /// <item><description>
    /// The distribution behaves similarly to the "RiskCumul" function in the Palisade's @Risk software.
    /// <see href="http://kb.palisade.com/index.php?pg=kb.page&id=51"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class EmpiricalDistribution : UnivariateDistributionBase, IBootstrappable
    {
       
        /// <summary>
        /// Constructs a Univariate Empirical CDF with default X {-0.5, 0, 0.5} and P values {0.1, 0.5, 0.9} and min = -1 and max = 1.
        /// </summary>
        public EmpiricalDistribution()
        {
            SetParameters(new[] { -0.5d, 0d, 0.5d }, new[] { 0.1d, 0.5d, 0.9d });
        }

        /// <summary>
        /// Constructs a Univariate Empirical CDF with specified parameters.
        /// </summary>
        /// <param name="XValues">Array of X values.</param>
        /// <param name="PValues">Array of probability values. Range 0 ≤ p ≤ 1.</param>
        public EmpiricalDistribution(IList<double> XValues, IList<double> PValues)
        {
            SetParameters(XValues, PValues);
          
        }

        /// <summary>
        /// Constructs a Univariate Empirical CDF with specified parameters.
        /// </summary>
        /// <param name="XValues">Array of X values.</param>
        /// <param name="PValues">Array of probability values. Range 0 ≤ p ≤ 1.</param>
        /// <param name="XOrder">Sort order of X values.</param>
        /// <param name="probabilityOrder">Sort order of probability values.</param>
        public EmpiricalDistribution(IList<double> XValues, IList<double> PValues, SortOrder XOrder, SortOrder probabilityOrder)
        {
            SetParameters(XValues, PValues, XOrder, probabilityOrder);

        }

        /// <summary>
        /// Constructs a Univariate Empirical CDF from ordered paired data.
        /// </summary>
        /// <param name="orderedPairedData">The ordered paired data.</param>
        public EmpiricalDistribution(OrderedPairedData orderedPairedData)
        {
            if (orderedPairedData.OrderX != SortOrder.Ascending) throw new ArgumentException("The x values must be in ascending order", nameof(orderedPairedData));
            opd = orderedPairedData;
            _xValues = orderedPairedData.Select(v => v.X).ToArray();

            // Probability values can be in ascending or descending order. 
            // If the sort order is "None", then we cannot use smart search.

            // Check if the probability values are actually in ascending order
            if (opd.OrderY == SortOrder.None)
            {            
                bool isAsc = true;
                _pValues = new double[opd.Count];
                for (int i = 0; i < opd.Count; i++)
                {
                    _pValues[i] = opd[i].Y;
                    if (i > 0 && _pValues[i] <= _pValues[i - 1])
                    {
                        isAsc = false;
                    }
                }
                if (isAsc == true)
                {
                    opd = new OrderedPairedData(_xValues, _pValues, true, SortOrder.Ascending, true, SortOrder.Ascending);
                }
                else
                {
                    opd.UseSmartSearch = false;
                }
            }
            else
            {
                _pValues = orderedPairedData.Select(v => v.Y).ToArray();
            }
        }

        /// <summary>
        /// Constructs a Univariate Empirical CDF from sample data.
        /// </summary>
        /// <param name="sample"></param>
        /// <param name="plottingPostionType">The plotting position formula type. Default = Weibull.</param>
        public EmpiricalDistribution(IList<double> sample, PlottingPositions.PlottingPostionType plottingPostionType = PlottingPositions.PlottingPostionType.Weibull)
        {
            _xValues = sample.ToArray();
            Array.Sort(_xValues);
            _pValues = PlottingPositions.Function(_xValues.Count(), plottingPostionType);
            opd = new OrderedPairedData(_xValues, _pValues, true, SortOrder.Ascending, true, SortOrder.Ascending);
        }

        private double[] _xValues;
        private double[] _pValues;
        private bool _parametersValid = true;
        private OrderedPairedData opd;
        private bool _momentsComputed = false;
        private double u1, u2, u3, u4;

        /// <summary>
        /// Returns the array of X values. Points On the cumulative curve are specified
        /// with increasing value and increasing probability.
        /// </summary>
        public ReadOnlyCollection<double> XValues => new ReadOnlyCollection<double>(_xValues);

        /// <summary>
        /// Returns the array of probability values. Points on the cumulative curve are specified
        /// with increasing value and increasing probability.
        /// </summary>
        public ReadOnlyCollection<double> ProbabilityValues => new ReadOnlyCollection<double>(_pValues);

        /// <summary>
        /// Returns the sort order of the X-values.
        /// </summary>
        public SortOrder XValueOrder => opd.OrderX;

        /// <summary>
        /// Returns the sort order of the Probability-values.
        /// </summary>
        public SortOrder ProbabilityOrder => opd.OrderY;

        /// <summary>
        /// Determines the interpolation transform for the X-values.
        /// </summary>
        public Transform XTransform { get; set; } = Transform.None;

        /// <summary>
        /// Determines the interpolation transform for the Probability-values.
        /// </summary>
        public Transform ProbabilityTransform { get; set; } = Transform.NormalZ;

        /// <summary>
        /// Returns the number of distribution parameters.
        /// </summary>
        public override int NumberOfParameters
        {
            get { return 2; }
        }

        /// <summary>
        /// Returns the continuous distribution type.
        /// </summary>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.Empirical; }
        }

        /// <summary>
        /// Returns the name of the distribution type as a string.
        /// </summary>
        public override string DisplayName
        {
            get { return "Univariate Empirical"; }
        }

        /// <summary>
        /// Returns the short display name of the distribution as a string.
        /// </summary>
        public override string ShortDisplayName
        {
            get { return "Uni. Emp"; }
        }

        /// <summary>
        /// Get distribution parameters in 2-column array of string.
        /// </summary>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                string Xstring = "{";
                string Pstring = "{";
                for (int i = 1; i < XValues.Count - 1; i++)
                {
                    Xstring += XValues[i].ToString();
                    Pstring += ProbabilityValues[i].ToString();
                    if (i < XValues.Count - 2)
                    {
                        Xstring += ",";
                        Pstring += ",";
                    }
                }
                Xstring += "}";
                Pstring += "}";
                parmString[0, 0] = "X Values";
                parmString[1, 0] = "P Values";
                parmString[0, 1] = Xstring;
                parmString[1, 1] = Pstring;
                return parmString;
            }
        }

        /// <summary>
        /// Gets the short form parameter names.
        /// </summary>
        public override string[] ParameterNamesShortForm
        {
            get { return new[] { "X()", "P()" }; }
        }

        /// <summary>
        /// Gets the full parameter names.
        /// </summary>
        public override string[] GetParameterPropertyNames
        {
            get { return new[] { nameof(XValues), nameof(ProbabilityValues) }; }
        }

        /// <summary>
        /// Get an array of parameters.
        /// </summary>
        public override double[] GetParameters
        {
            get { throw new NotImplementedException(); }
        }

        /// <summary>
        /// Determines whether the parameters are valid or not.
        /// </summary>
        public override bool ParametersValid
        {
            get { return _parametersValid; }
        }

        /// <summary>
        /// Compute central moments of the distribution.
        /// </summary>
        private void ComputeMoments()
        {
            var mom = CentralMoments(200);
            u1 = mom[0];
            u2 = mom[1];
            u3 = mom[2];
            u4 = mom[3];
            _momentsComputed = true;
        }

        /// <summary>
        /// Gets the mean of the distribution.
        /// </summary>
        public override double Mean
        {
            get
            {
                if (!_momentsComputed) ComputeMoments();
                return u1;
            }
        }

        /// <summary>
        /// Gets the median of the distribution.
        /// </summary>
        public override double Median
        {
            get { return InverseCDF(0.5d); }
        }

        /// <summary>
        /// Gets the mode of the distribution.
        /// </summary>
        public override double Mode
        {
            get { return double.NaN; }
        }

        /// <summary>
        /// Gets the standard deviation of the distribution.
        /// </summary>
        public override double StandardDeviation
        {
            get
            {
                if (!_momentsComputed) ComputeMoments();
                return u2;
            }
        }

        /// <summary>
        /// Gets the skew of the distribution.
        /// </summary>
        public override double Skew
        {
            get
            {
                if (!_momentsComputed) ComputeMoments();
                return u3;
            }
        }

        /// <summary>
        /// Gets the kurtosis of the distribution.
        /// </summary>
        public override double Kurtosis
        {
            get
            {
                if (!_momentsComputed) ComputeMoments();
                return u4;
            }
        }

        /// <summary>
        /// Gets the minimum of the distribution.
        /// </summary>
        public override double Minimum
        {
            get { return XValues.First(); }
        }

        /// <summary>
        /// Gets the maximum of the distribution.
        /// </summary>
        public override double Maximum
        {
            get { return XValues.Last(); }
        }

        /// <summary>
        /// Gets the minimum values allowable for each parameter.
        /// </summary>
        public override double[] MinimumOfParameters
        {
            get { return new double[] { double.MinValue, 0d }; }
        }

        /// <summary>
        /// Gets the maximum values allowable for each parameter.
        /// </summary>
        public override double[] MaximumOfParameters
        {
            get { return new double[] { double.MaxValue, 1d }; }
        }

        /// <summary>
        /// Bootstrap the distribution based on a sample size and parameter estimation method.
        /// </summary>
        /// <param name="estimationMethod">The parameter estimation method.</param>
        /// <param name="sampleSize">Size of the random sample to generate.</param>
        /// <param name="seed">Optional. Seed for random number generator. Default = 12345.</param>
        /// <returns>
        /// Returns a bootstrapped distribution.
        /// </returns>
        public IUnivariateDistribution Bootstrap(ParameterEstimationMethod estimationMethod, int sampleSize, int seed = 12345)
        {
            var sample = GenerateRandomValues(seed, sampleSize);
            return new EmpiricalDistribution(sample);
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="xValues">Array of X values.</param>
        /// <param name="pValues">Array of probability values. Range 0 ≤ p ≤ 1.</param>
        public void SetParameters(IList<double> xValues, IList<double> pValues)
        {
            _xValues = xValues.ToArray();
            _pValues = pValues.ToArray();
            opd = new OrderedPairedData(xValues, pValues, true, SortOrder.Ascending, true, SortOrder.Ascending);
            _momentsComputed = false;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="xValues">Array of X values.</param>
        /// <param name="pValues">Array of probability values. Range 0 ≤ p ≤ 1.</param>
        /// <param name="XOrder">Sort order of X values.</param>
        /// <param name="probabilityOrder">Sort order of probability values.</param>
        public void SetParameters(IList<double> xValues, IList<double> pValues, SortOrder XOrder, SortOrder probabilityOrder)
        {
            _xValues = xValues.ToArray();
            _pValues = pValues.ToArray();
            opd = new OrderedPairedData(xValues, pValues, true, XOrder, true, probabilityOrder);
            _momentsComputed = false;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="parameters">A list of parameters.</param>
        public override void SetParameters(IList<double> parameters)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="parameters">A list of parameters.</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            return null;
        }

        /// <summary>
        /// Gets the Probability Density Function (PDF) of the distribution evaluated at a point X.
        /// </summary>
        /// <param name="X">A single point in the distribution range.</param>
        public override double PDF(double X)
        {
            if (X < Minimum || X > Maximum) return 0.0d;
            var d = NumericalDerivative.Derivative(CDF, X);
            return d < 0d ? 0d : d;
        }

        /// <summary>
        /// Returns the Probability Density Function (PDF) of the distribution.
        /// </summary>
        /// <param name="xl">Lower x value.</param>
        /// <param name="xu">Upper x value</param>
        /// <returns></returns>
        public double PDF(double xl, double xu)
        {
            if (xu == xl) return PDF(xu);
            if (xu < xl) return 0.0d;
            if (xl < Minimum || xl > Maximum || xu < Minimum || xu > Maximum) return 0.0d;
            var d = (CDF(xu) - CDF(xl)) / (xu - xl);
            return d < 0d ? 0d : d;
        }

        /// <summary>
        /// Gets the Cumulative Distribution Function (CDF) for the distribution evaluated at a point X.
        /// </summary>
        /// <param name="X">A single point in the distribution range.</param>
        public override double CDF(double X)
        {
            //if (X < Minimum) return 0d;
            //if (X > Maximum) return 1d;
            double p = 0;
            if (opd.OrderY == SortOrder.Ascending || opd.OrderY == SortOrder.None)
            {
                p = opd.Interpolate(X, true, XTransform, ProbabilityTransform);
            }
            else
            {
                // If descending then it is a survival function
                p = 1d - opd.Interpolate(X, true, XTransform, ProbabilityTransform);
            }
            return p < 0d ? 0d : p > 1d ? 1d : p;
        }

        /// <summary>
        /// Gets the Inverse Cumulative Distribution Function (ICFD) of the distribution evaluated at a probability.
        /// </summary>
        /// <param name="probability">Probability between 0 and 1.</param>
        /// <returns>
        /// Returns for a given probability in the probability distribution of a random variable,
        /// the value at which the probability of the random variable is less than or equal to the
        /// given probability.
        /// </returns>
        /// <remarks>
        /// This function is also know as the Quantile Function.
        /// </remarks>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (probability < 0.0d || probability > 1.0d)
                throw new ArgumentOutOfRangeException(nameof(probability), "Probability must be between 0 and 1.");

            double min = Minimum;
            double max = Maximum;
            if (probability <= 1E-16) return min;
            if (probability >= 1 - 1E-16) return max;
            double x = 0;
            if (opd.OrderY == SortOrder.Ascending || opd.OrderY == SortOrder.None)
            {
                x = opd.Interpolate(probability, false, XTransform, ProbabilityTransform);
            }
            else
            {
                // If descending then it is a survival function
                x = opd.Interpolate(1d - probability, false, XTransform, ProbabilityTransform);
            }
            return x < min ? min : x > max ? max : x;
        }
       
        /// <summary>
        /// Creates a copy of the distribution.
        /// </summary>
        public override UnivariateDistributionBase Clone()
        {
            return new EmpiricalDistribution(XValues, ProbabilityValues) { XTransform = XTransform, ProbabilityTransform = ProbabilityTransform };
        }

    }
}