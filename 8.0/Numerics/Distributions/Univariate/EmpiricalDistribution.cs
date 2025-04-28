/*
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using Numerics.Data;
using Numerics.Data.Statistics;
using Numerics.Mathematics;
using Numerics.Mathematics.Optimization;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Univariate Empirical distribution.
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
    /// This distribution specifies a cumulative distribution with n points. The range of the distribution
    /// is set by the minimum and maximum arguments. Each point on the cumulative curve has an X value and
    /// a probability. Points on the cumulative curve must be entered with increasing value and increasing
    /// probability. Even though the (X,p) points define the distribution, and value between the minimum
    /// and maximum can be returned.
    /// </para>
    /// <para>
    /// <b> References:</b>
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
            SetParameters([-0.5d, 0d, 0.5d], [0.1d, 0.5d, 0.9d]);
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
            _momentsComputed = false;
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
            _momentsComputed = false;
        }

        private double[] _xValues;
        private double[] _pValues;
        private OrderedPairedData opd;
        private bool _momentsComputed = false;
        private double u1, u2, u3, u4;

        /// <summary>
        /// Returns the array of X values. Points On the cumulative curve are specified
        /// with increasing value and increasing probability.
        /// </summary>
        public ReadOnlyCollection<double> XValues => new(_xValues);

        /// <summary>
        /// Returns the array of probability values. Points on the cumulative curve are specified
        /// with increasing value and increasing probability.
        /// </summary>
        public ReadOnlyCollection<double> ProbabilityValues => new(_pValues);

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

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 2; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.Empirical; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Univariate Empirical"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "Uni. Emp"; }
        }

        /// <inheritdoc/>
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

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["X()", "P()"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(XValues), nameof(ProbabilityValues)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { throw new NotImplementedException(); }
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

        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                if (!_momentsComputed) ComputeMoments();
                return u1;
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
            get
            {
                var brent = new BrentSearch(PDF, InverseCDF(0.001), InverseCDF(0.999));
                brent.Maximize();
                return brent.BestParameterSet.Values[0];
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                if (!_momentsComputed) ComputeMoments();
                return u2;
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {
                if (!_momentsComputed) ComputeMoments();
                return u3;
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                if (!_momentsComputed) ComputeMoments();
                return u4;
            }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return XValues.First(); }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return XValues.Last(); }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [double.MinValue, 0d]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.MaxValue, 1d]; }
        }

        /// <inheritdoc/>
        public IUnivariateDistribution Bootstrap(ParameterEstimationMethod estimationMethod, int sampleSize, int seed = -1)
        {
            var sample = GenerateRandomValues(sampleSize, seed);
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

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            throw new NotImplementedException();
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            return null;
        }

        /// <inheritdoc/>
        public override double PDF(double X)
        {
            if (X < Minimum || X > Maximum) return 0.0d;

            double dFdx = 0;
            double x = X;
            double h = NumericalDerivative.CalculateStepSize(X);
            if (X <= _xValues.First())
            {
                // If x is on the boundary, use a one-sided two-point difference
                x = _xValues.First();
                h = NumericalDerivative.CalculateStepSize(x);
                dFdx = (CDF(x + h) - CDF(x)) / h;
            }
            else if (X >= _xValues.Last())
            {
                // If x is on the boundary, use a one-sided two-point difference
                x = _xValues.Last();
                h = NumericalDerivative.CalculateStepSize(x);
                dFdx = (CDF(x) - CDF(x - h)) / h;
            }
            else
            {
                // otherwise central difference
                dFdx = (CDF(x + h) - CDF(x - h)) / (2 * h);
            }
            return dFdx < 0 ? 0 : dFdx;
        }

        /// <summary>
        /// Returns the Probability Density Function (PDF) of the distribution.
        /// </summary>
        /// <param name="xl">Lower x value.</param>
        /// <param name="xu">Upper x value</param>
        public double PDF(double xl, double xu)
        {
            if (xu == xl) return PDF(xu);
            if (xu < xl) return 0.0d;
            if (xl < Minimum || xl > Maximum || xu < Minimum || xu > Maximum) return 0.0d;
            var d = (CDF(xu) - CDF(xl)) / (xu - xl);
            return d < 0d ? 0d : d;
        }

        /// <inheritdoc/>
        public override double CDF(double X)
        {
            double p = 0;
            if (opd.OrderY == SortOrder.Ascending || opd.OrderY == SortOrder.None)
            {
                p = opd.GetYFromX(X, XTransform, ProbabilityTransform);
            }
            else
            {
                // If descending then it is a survival function
                p = 1d - opd.GetYFromX(X, XTransform, ProbabilityTransform);
            }
            return p < 0d ? 0d : p > 1d ? 1d : p;
        }

        /// <inheritdoc/>
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
                x = opd.GetXFromY(probability, XTransform, ProbabilityTransform);
            }
            else
            {
                // If descending then it is a survival function
                x = opd.GetXFromY(1d - probability, XTransform, ProbabilityTransform);
            }
            return x < min ? min : x > max ? max : x;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new EmpiricalDistribution(XValues, ProbabilityValues) { XTransform = XTransform, ProbabilityTransform = ProbabilityTransform };
        }

    }
}