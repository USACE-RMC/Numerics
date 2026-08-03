using Numerics.Data;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.RootFinding;
using Numerics.Sampling;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Xml.Linq;

namespace Numerics.Distributions
{
    /// <summary>
    /// A Mixture distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    /// When zero inflation is enabled, the distribution is a positive-hurdle mixture:
    /// the zero atom has probability <see cref="ZeroWeight"/>, component weights sum to
    /// the remaining mass, and every component is conditioned on a value strictly greater
    /// than zero. <see cref="PDF(double)"/> at zero reports the atom probability under the
    /// mixed Lebesgue-plus-Dirac reference measure.
    /// </para>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Mixture_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class Mixture : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IBootstrappable
    {

        /// <summary>
        /// Construct new mixture distribution.
        /// </summary>
        /// <param name="weights">The mixture weights.</param>
        /// <param name="distributions">The mixture distributions.</param>
        public Mixture(double[] weights, UnivariateDistributionBase[] distributions)
        {
            SetParameters(weights, distributions);
        }

        /// <summary>
        /// Construct new mixture distribution.
        /// </summary>
        /// <param name="weights">The mixture weights.</param>
        /// <param name="distributions">The mixture distributions.</param>
        public Mixture(double[] weights, IUnivariateDistribution[] distributions)
        {
            SetParameters(weights, distributions);
        }

        private double[] _weights = null!;
        private bool _isZeroInflated;
        private double _zeroWeight;
        private UnivariateDistributionBase[] _distributions = null!;
        private EmpiricalDistribution _empiricalCDF = null!;
        private bool _momentsComputed = false;
        private double u1, u2, u3, u4;
        private bool _empiricalCDFCreated = false;

        /// <summary>
        /// Returns the array of distribution weights.
        /// </summary>
        public double[] Weights => _weights;

        /// <summary>
        /// Returns the array of univariate probability distributions.
        /// </summary>
        public UnivariateDistributionBase[] Distributions => _distributions;

        /// <summary>
        /// Gets or sets whether a separate probability weight is assigned to values less than or equal to zero.
        /// </summary>
        /// <remarks>
        /// Enabling zero inflation rescales finite, nonnegative component weights so their sum
        /// equals <c>1 - <see cref="ZeroWeight"/></c>.
        /// </remarks>
        public bool IsZeroInflated
        {
            get { return _isZeroInflated; }
            set
            {
                _isZeroInflated = value;
                if (_isZeroInflated)
                {
                    NormalizeComponentWeights();
                }
                RefreshConfigurationState();
            }
        }

        /// <summary>
        /// Gets or sets the zero-value probability weight used when the mixture is zero-inflated.
        /// </summary>
        /// <remarks>
        /// When zero inflation is enabled, assigning this property rescales finite, nonnegative
        /// component weights to the remaining probability mass.
        /// </remarks>
        public double ZeroWeight
        {
            get { return _zeroWeight; }
            set
            {
                _zeroWeight = value;
                if (IsZeroInflated)
                {
                    NormalizeComponentWeights();
                }
                RefreshConfigurationState();
            }
        }

        /// <summary>
        /// Rescales valid component weights to the probability mass remaining after zero inflation.
        /// </summary>
        /// <remarks>
        /// Invalid zero weights and invalid component weights are left unchanged so parameter
        /// validation can report the original configuration error.
        /// </remarks>
        private void NormalizeComponentWeights()
        {
            if (_weights is null || _weights.Length == 0 ||
                double.IsNaN(ZeroWeight) || double.IsInfinity(ZeroWeight) ||
                ZeroWeight < 0.0 || ZeroWeight > 1.0)
            {
                return;
            }

            double sum = 0.0;
            for (int i = 0; i < _weights.Length; i++)
            {
                if (double.IsNaN(_weights[i]) || double.IsInfinity(_weights[i]) || _weights[i] < 0.0)
                {
                    return;
                }
                sum += _weights[i];
            }

            if (sum <= 0.0 || double.IsInfinity(sum))
            {
                return;
            }

            double scale = (1.0 - ZeroWeight) / sum;
            for (int i = 0; i < _weights.Length; i++)
            {
                _weights[i] *= scale;
            }
        }

        /// <summary>
        /// Determines whether a value is finite on every target framework.
        /// </summary>
        /// <param name="value">The value to inspect.</param>
        /// <returns><see langword="true"/> when the value is neither NaN nor infinite.</returns>
        private static bool IsFinite(double value)
        {
            return !double.IsNaN(value) && !double.IsInfinity(value);
        }

        /// <summary>
        /// Restricts a value to an inclusive interval on every target framework.
        /// </summary>
        /// <param name="value">The value to restrict.</param>
        /// <param name="minimum">The inclusive lower bound.</param>
        /// <param name="maximum">The inclusive upper bound.</param>
        /// <returns>The restricted value.</returns>
        private static double Clamp(double value, double minimum, double maximum)
        {
            return value < minimum ? minimum : value > maximum ? maximum : value;
        }

        /// <summary>
        /// Returns the next representable value greater than the supplied value.
        /// </summary>
        /// <param name="value">The starting value.</param>
        /// <returns>The adjacent representable value toward positive infinity.</returns>
        private static double BitIncrement(double value)
        {
            if (double.IsNaN(value) || value == double.PositiveInfinity) return value;
            if (value == 0.0) return double.Epsilon;
            long bits = BitConverter.DoubleToInt64Bits(value);
            return BitConverter.Int64BitsToDouble(value > 0.0 ? bits + 1 : bits - 1);
        }

        /// <summary>
        /// Returns the next representable value less than the supplied value.
        /// </summary>
        /// <param name="value">The starting value.</param>
        /// <returns>The adjacent representable value toward negative infinity.</returns>
        private static double BitDecrement(double value)
        {
            if (double.IsNaN(value) || value == double.NegativeInfinity) return value;
            if (value == 0.0) return -double.Epsilon;
            long bits = BitConverter.DoubleToInt64Bits(value);
            return BitConverter.Int64BitsToDouble(value > 0.0 ? bits - 1 : bits + 1);
        }

        /// <summary>
        /// Converts a unit-interval draw to a component CDF probability conditional on a positive value.
        /// </summary>
        /// <param name="componentIndex">The zero-based component index.</param>
        /// <param name="conditionalProbability">The probability on the positive-conditional scale.</param>
        /// <returns>A component CDF probability strictly above the CDF at zero and strictly below one.</returns>
        private double PositiveConditionalQuantileProbability(int componentIndex, double conditionalProbability)
        {
            TryGetPositiveMass(componentIndex, out double positiveMass);
            double cdfAtZero = Distributions[componentIndex].CDF(0.0);
            double probability = cdfAtZero + conditionalProbability * positiveMass;
            if (probability <= cdfAtZero) probability = BitIncrement(cdfAtZero);
            if (probability >= 1.0) probability = BitDecrement(1.0);
            return probability;
        }
        /// <summary>
        /// Gets the probability that a component produces a strictly positive value.
        /// </summary>
        /// <param name="componentIndex">The zero-based component index.</param>
        /// <param name="positiveMass">The strictly positive probability mass.</param>
        /// <returns><see langword="true"/> when the mass is finite and positive; otherwise, <see langword="false"/>.</returns>
        private bool TryGetPositiveMass(int componentIndex, out double positiveMass)
        {
            positiveMass = Distributions[componentIndex].CCDF(0.0);
            return IsFinite(positiveMass) && positiveMass > 0.0;
        }

        /// <summary>
        /// Evaluates a component density conditional on a strictly positive value.
        /// </summary>
        /// <param name="componentIndex">The zero-based component index.</param>
        /// <param name="x">The value at which to evaluate the density.</param>
        /// <returns>The positive-conditional density.</returns>
        private double PositiveConditionalPDF(int componentIndex, double x)
        {
            return x > 0.0 && TryGetPositiveMass(componentIndex, out double positiveMass)
                ? Distributions[componentIndex].PDF(x) / positiveMass
                : 0.0;
        }

        /// <summary>
        /// Evaluates a component log density conditional on a strictly positive value.
        /// </summary>
        /// <param name="componentIndex">The zero-based component index.</param>
        /// <param name="x">The value at which to evaluate the log density.</param>
        /// <returns>The positive-conditional log density.</returns>
        private double PositiveConditionalLogPDF(int componentIndex, double x)
        {
            return x > 0.0 && TryGetPositiveMass(componentIndex, out double positiveMass)
                ? Distributions[componentIndex].LogPDF(x) - Math.Log(positiveMass)
                : double.NegativeInfinity;
        }

        /// <summary>
        /// Evaluates a component distribution function conditional on a strictly positive value.
        /// </summary>
        /// <param name="componentIndex">The zero-based component index.</param>
        /// <param name="x">The value at which to evaluate the distribution function.</param>
        /// <returns>The positive-conditional cumulative probability.</returns>
        private double PositiveConditionalCDF(int componentIndex, double x)
        {
            if (x <= 0.0 || !TryGetPositiveMass(componentIndex, out double positiveMass)) return 0.0;
            double probability = (Distributions[componentIndex].CDF(x) - Distributions[componentIndex].CDF(0.0)) / positiveMass;
            return Clamp(probability, 0.0, 1.0);
        }

        /// <summary>
        /// Evaluates a component survival function conditional on a strictly positive value.
        /// </summary>
        /// <param name="componentIndex">The zero-based component index.</param>
        /// <param name="x">The value at which to evaluate the survival function.</param>
        /// <returns>The positive-conditional survival probability.</returns>
        private double PositiveConditionalCCDF(int componentIndex, double x)
        {
            if (x < 0.0) return 1.0;
            if (!TryGetPositiveMass(componentIndex, out double positiveMass)) return double.NaN;
            double probability = Distributions[componentIndex].CCDF(x) / positiveMass;
            return Clamp(probability, 0.0, 1.0);
        }
        /// <summary>
        /// Refreshes validity and cached results after zero-inflation configuration changes.
        /// </summary>
        private void RefreshConfigurationState()
        {
            if (_weights is null || _distributions is null)
            {
                _parametersValid = false;
            }
            else
            {
                _parametersValid = ValidateParameters(GetParameters, false) is null;
            }
            _momentsComputed = false;
            _empiricalCDFCreated = false;
        }

        /// <summary>
        /// Determines the interpolation transform for the X-values.
        /// </summary>
        public Transform XTransform { get; set; } = Transform.None;

        /// <summary>
        /// Determines the interpolation transform for the Probability-values.
        /// </summary>
        public Transform ProbabilityTransform { get; set; } = Transform.NormalZ;

        /// <summary>
        /// The maximum iterations in the Expectation Maximization algorithm. Default = 1,000. 
        /// </summary>
        public int MaxIterations { get; set; } = 1000;

        /// <summary>
        /// The relative tolerance for convergence. Default = 1E-8.
        /// </summary>
        public double Tolerance { get; set; } = 1E-8;

        /// <summary>
        /// The total number of iterations required to find the MLE.
        /// </summary>
        public int Iterations { get; private set; }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get
            {
                int sum = 0;
                sum += Distributions.Count();
                for (int i = 0; i < Distributions.Count(); i++)
                    sum += Distributions[i].NumberOfParameters;
                return sum;
            }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type => UnivariateDistributionType.Mixture;

        /// <inheritdoc/>
        public override string DisplayName => "Mixture";

        /// <inheritdoc/>
        public override string ShortDisplayName => "MIX";

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                string Wstring = "{";
                string Dstring = "{";
                for (int i = 0; i < Weights.Count(); i++)
                {
                    Wstring += Weights[i].ToString();
                    Dstring += Distributions[i].DisplayName;
                    if (i < Weights.Count() - 1)
                    {
                        Wstring += ",";
                        Dstring += ",";
                    }
                }
                Wstring += "}";
                Dstring += "}";
                parmString[0, 0] = "Weights";
                parmString[1, 0] = "Distributions";
                parmString[0, 1] = Wstring;
                parmString[1, 1] = Dstring;
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNames
        {
            get
            {
                var result = new List<string>();
                for (int i = 1; i <= Distributions.Count(); i++)
                {
                    result.Add("Weight " + i.ToString());
                }
                for (int i = 0; i < Distributions.Count(); i++)
                {
                    for (int j = 0; j < Distributions[i].ParameterNames.Length; j++)
                    {
                        result.Add("D" + (i + 1).ToString() + " " + Distributions[i].ParameterNames[j]);
                    }
                }
                return result.ToArray();
            }
        }


        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get
            {
                var result = new List<string>();
                for (int i = 1; i <= Distributions.Count(); i++)
                {
                    result.Add("W" + i.ToString());
                }
                for (int i = 0; i < Distributions.Count(); i++)
                {
                    for (int j = 0; j < Distributions[i].ParameterNamesShortForm.Length; j++)
                    {
                        result.Add("D" + (i + 1).ToString() + " " + Distributions[i].ParameterNamesShortForm[j]);
                    }
                }
                return result.ToArray();
            }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get
            {
                var result = new List<double>();
                result.AddRange(Weights);
                for (int i = 0; i < Distributions.Count(); i++)
                {
                    result.AddRange(Distributions[i].GetParameters);
                }                  
                return result.ToArray();
            }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Weights), nameof(Distributions)]; }
        }

        /// <summary>
        /// Compute central moments of the distribution.
        /// </summary>
        private void ComputeMoments()
        {
            var mom = CentralMoments(1000);
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
                if (!_momentsComputed) 
                    ComputeMoments();
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
                if (!_momentsComputed) 
                    ComputeMoments();
                return u2;
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {
                if (!_momentsComputed) 
                    ComputeMoments();
                return u3;
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                if (!_momentsComputed) 
                    ComputeMoments();
                return u4;
            }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return IsZeroInflated ? 0.0 : Distributions.Min(p => p.Minimum); }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return Distributions.Max(p => p.Maximum); }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        { 
            get 
            {
                var result = new List<double>();
                if (IsZeroInflated) { result.Add(0.0); }
                for (int i = 0; i < Distributions.Count(); i++)
                {
                    result.Add(0.0);
                }
                for (int i = 0; i < Distributions.Count(); i++)
                {
                    result.AddRange(Distributions[i].MinimumOfParameters);
                }
                return result.ToArray(); 
            }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get
            {
                var result = new List<double>();
                if (IsZeroInflated) { result.Add(1.0); }
                for (int i = 0; i < Distributions.Count(); i++)
                {
                    result.Add(1.0);
                }
                for (int i = 0; i < Distributions.Count(); i++)
                {
                    result.AddRange(Distributions[i].MaximumOfParameters);
                }
                return result.ToArray();
            }
        }

        /// <inheritdoc/>
        public void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod)
        {
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
            var newDistribution = (Mixture)Clone();
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="weights">The mixture weights.</param>
        /// <param name="distributions">The mixture distributions.</param>
        public void SetParameters(double[] weights, UnivariateDistributionBase[] distributions)
        {
            if (weights == null) throw new ArgumentNullException(nameof(Weights));
            if (distributions == null) throw new ArgumentNullException(nameof(Distributions));
            if (weights.Length != distributions.Length)
                throw new ArgumentException("The weight and distribution arrays must have the same length.", nameof(Weights));

            _weights = weights.ToArray();
            _distributions = distributions.ToArray();
            _parametersValid = ValidateParameters(GetParameters, false) is null;
            _momentsComputed = false;
            _empiricalCDFCreated = false;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="weights">The mixture weights.</param>
        /// <param name="distributions">The mixture distributions.</param>
        public void SetParameters(double[] weights, IUnivariateDistribution[] distributions)
        {
            if (weights == null) throw new ArgumentNullException(nameof(Weights));
            if (distributions == null) throw new ArgumentNullException(nameof(Distributions));
            if (weights.Length != distributions.Length)
                throw new ArgumentException("The weight and distribution arrays must have the same length.", nameof(Weights));

            _weights = weights.ToArray();
            _distributions = new UnivariateDistributionBase[distributions.Length];
            for (int i = 0; i < distributions.Length; i++)
            {
                _distributions[i] = (UnivariateDistributionBase)distributions[i];
            }
            _parametersValid = ValidateParameters(GetParameters, false) is null;
            _momentsComputed = false;
            _empiricalCDFCreated = false;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="weights">The mixture weights.</param>
        /// <param name="parameters">The mixture distribution parameters.</param>
        public void SetParameters(double[] weights, double[] parameters)
        {
            if (weights == null) throw new ArgumentNullException(nameof(Weights));
            if (parameters == null) throw new ArgumentNullException(nameof(parameters));
            if (weights.Length != Distributions.Length)
                throw new ArgumentException("The weight and distribution arrays must have the same length.", nameof(Weights));
            if (parameters.Length != Distributions.Sum(x => x.NumberOfParameters))
            {
                throw new ArgumentException("The length of the parameter array is invalid.", nameof(parameters));
            }

            double[] parameterCopy = parameters.ToArray();

            // Set weights
            _weights = weights.ToArray();
            // Set distribution parameters
            int t = 0;
            for (int i = 0; i < Distributions.Count(); i++)
            {
                var parms = new List<double>();
                for (int j = t; j < t + Distributions[i].NumberOfParameters; j++)
                {
                    parms.Add(parameterCopy[j]);
                }
                Distributions[i].SetParameters(parms);
                t += Distributions[i].NumberOfParameters;
            }
            _parametersValid = ValidateParameters(GetParameters, false) is null;
            _momentsComputed = false;
            _empiricalCDFCreated = false;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            if (parameters == null) throw new ArgumentNullException(nameof(parameters));
            if (parameters.Count != NumberOfParameters)
            {
                throw new ArgumentException("The length of the parameter array is invalid.", nameof(parameters));
            }

            double[] parameterCopy = parameters.ToArray();

            // Set the weights.
            int parameterIndex = 0;
            for (int i = 0; i < Distributions.Count(); i++)
            {
                Weights[i] = parameterCopy[parameterIndex++];
            }

            // Set the distribution parameters.
            for (int i = 0; i < Distributions.Count(); i++)
            {
                double[] distributionParameters = parameterCopy
                    .Skip(parameterIndex)
                    .Take(Distributions[i].NumberOfParameters)
                    .ToArray();
                Distributions[i].SetParameters(distributionParameters);
                parameterIndex += Distributions[i].NumberOfParameters;
            }

            _parametersValid = ValidateParameters(GetParameters, false) is null;
            _momentsComputed = false;
            _empiricalCDFCreated = false;
        }

        /// <summary>
        /// Set the distribution parameters from a referenced array. Weights are normalized to the configured simplex.
        /// </summary>
        /// <param name="parameters">The array of parameters. The caller's array is not modified.</param>
        public void SetParameters(ref double[] parameters)
        {
            if (parameters == null) return;
            if (Weights == null || Weights.Length == 0) return;
            if (Distributions == null || Distributions.Count() == 0) return;

            double[] parameterCopy = parameters.ToArray();
            if (Distributions.Count() == 1 && parameterCopy.Length == Distributions[0].NumberOfParameters)
            {
                Weights[0] = IsZeroInflated ? 1.0 - ZeroWeight : 1.0;
                Distributions[0].SetParameters(parameterCopy);
            }
            else
            {
                int componentCount = Distributions.Count();
                int parameterIndex = componentCount;
                double weightSum = 0.0;

                for (int i = 0; i < componentCount; i++)
                {
                    Weights[i] = parameterCopy[i];
                    weightSum += Weights[i];
                }

                double componentMass = IsZeroInflated ? 1.0 - ZeroWeight : 1.0;
                if (weightSum <= 0.0 || !IsFinite(weightSum))
                {
                    double uniformWeight = componentMass / componentCount;
                    for (int i = 0; i < componentCount; i++) Weights[i] = uniformWeight;
                }
                else
                {
                    double scale = componentMass / weightSum;
                    for (int i = 0; i < componentCount; i++) Weights[i] *= scale;
                }

                for (int i = 0; i < componentCount; i++)
                {
                    double[] distributionParameters = parameterCopy
                        .Skip(parameterIndex)
                        .Take(Distributions[i].NumberOfParameters)
                        .ToArray();
                    Distributions[i].SetParameters(distributionParameters);
                    parameterIndex += Distributions[i].NumberOfParameters;
                }
            }

            _parametersValid = ValidateParameters(GetParameters, false) is null;
            _momentsComputed = false;
            _empiricalCDFCreated = false;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (IsZeroInflated && (!IsFinite(ZeroWeight) || ZeroWeight < 0.0 || ZeroWeight >= 1.0))
            {
                var exception = new ArgumentOutOfRangeException(
                    nameof(ZeroWeight),
                    "The zero value weight must be finite and greater than or equal to 0 and less than 1.");
                if (throwException) throw exception;
                return exception;
            }

            for (int i = 0; i < Distributions.Count(); i++)
            {
                if (!IsFinite(Weights[i]) || Weights[i] < 0.0 || Weights[i] > 1.0)
                {
                    var exception = new ArgumentOutOfRangeException(
                        nameof(Weights),
                        "The weights must be finite and between 0 and 1.");
                    if (throwException) throw exception;
                    return exception;
                }
            }

            double totalMass = IsZeroInflated ? ZeroWeight : 0.0;
            for (int i = 0; i < Distributions.Count(); i++) totalMass += Weights[i];
            if (!IsFinite(totalMass) || !totalMass.AlmostEquals(1.0, 1E-8))
            {
                var exception = new ArgumentOutOfRangeException(
                    nameof(Weights),
                    IsZeroInflated
                        ? "The component weights must sum to 1 minus the zero value weight."
                        : "The weights must sum to 1.0.");
                if (throwException) throw exception;
                return exception;
            }

            for (int i = 0; i < Distributions.Count(); i++)
            {
                if (!Distributions[i].ParametersValid)
                {
                    var exception = new ArgumentOutOfRangeException(
                        nameof(Distributions),
                        "Distribution " + (i + 1).ToString() + " has invalid parameters.");
                    if (throwException) throw exception;
                    return exception;
                }

                if (IsZeroInflated && !TryGetPositiveMass(i, out _))
                {
                    var exception = new ArgumentOutOfRangeException(
                        nameof(Distributions),
                        "Distribution " + (i + 1).ToString() + " must have finite, positive probability above zero.");
                    if (throwException) throw exception;
                    return exception;
                }
            }

            return null;
        }

        /// <inheritdoc/>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            var initialVals = new double[NumberOfParameters];
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];

            // Weights are first
            int t = 0;
            for (int i = 0; i < Distributions.Count(); i++)
            {
                initialVals[i] = IsZeroInflated ? (1d - ZeroWeight) / Distributions.Count() : 1d / Distributions.Count();
                lowerVals[i] = 0.0;
                upperVals[i] = 1.0;
                t += 1;
            }

            for (int i = 0; i < Distributions.Count(); i++)
            {
                var tuple = ((IMaximumLikelihoodEstimation)Distributions[i]).GetParameterConstraints(sample);
                var initials = tuple.Item1;
                var lowers = tuple.Item2;
                var uppers = tuple.Item3;

                for (int j = t; j < t + Distributions[i].NumberOfParameters; j++)
                {
                    initialVals[j] = initials[j - t];
                    lowerVals[j] = lowers[j - t];
                    upperVals[j] = uppers[j - t];
                }
                t += Distributions[i].NumberOfParameters;
            }
            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
        }

        /// <inheritdoc/>
        public double[] MLE(IList<double> sample)
        {
            ValidateParameters(GetParameters, true);

            int observationCount = sample.Count;
            int distributionParameterCount = Distributions.Sum(x => x.NumberOfParameters);
            int componentCount = Distributions.Count();

            if (IsZeroInflated)
            {
                for (int rowIndex = 0; rowIndex < observationCount; rowIndex++)
                {
                    if (sample[rowIndex] < 0.0)
                    {
                        throw new InvalidOperationException(
                            "Mixture EM row " + rowIndex.ToString(CultureInfo.InvariantCulture) +
                            " has negative exact value " + sample[rowIndex].ToString("R", CultureInfo.InvariantCulture) +
                            " in a zero-inflated model.");
                    }
                }
            }

            Tuple<double[], double[], double[]> constraints = GetParameterConstraints(sample);
            double[] initialParameters = constraints.Item1.Subset(componentCount);
            double[] lowerParameters = constraints.Item2.Subset(componentCount);
            double[] upperParameters = constraints.Item3.Subset(componentCount);

            double[] mleWeights = constraints.Item1.Subset(0, componentCount - 1);
            double[] mleParameters = initialParameters;
            var responsibilities = new double[observationCount, componentCount];
            double oldLogLikelihood = double.MinValue;
            double newLogLikelihood = double.MinValue;

            double EStep(double[] parameters)
            {
                var distribution = (Mixture)Clone();
                distribution.SetParameters(mleWeights, parameters);
                double logLikelihood = 0.0;

                for (int rowIndex = 0; rowIndex < observationCount; rowIndex++)
                {
                    double value = sample[rowIndex];
                    if (IsZeroInflated && value == 0.0)
                    {
                        if (!IsFinite(ZeroWeight) || ZeroWeight <= 0.0)
                        {
                            throw CreateImpossibleRowException(rowIndex, value);
                        }

                        for (int componentIndex = 0; componentIndex < componentCount; componentIndex++)
                        {
                            responsibilities[rowIndex, componentIndex] = 0.0;
                        }
                        logLikelihood += Math.Log(ZeroWeight);
                        continue;
                    }

                    double maximumLogProbability = double.NegativeInfinity;
                    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++)
                    {
                        double componentLogDensity = IsZeroInflated
                            ? distribution.PositiveConditionalLogPDF(componentIndex, value)
                            : distribution.Distributions[componentIndex].LogPDF(value);
                        double logProbability = Math.Log(mleWeights[componentIndex]) + componentLogDensity;
                        responsibilities[rowIndex, componentIndex] = logProbability;
                        if (logProbability > maximumLogProbability) maximumLogProbability = logProbability;
                    }

                    if (!IsFinite(maximumLogProbability))
                    {
                        throw CreateImpossibleRowException(rowIndex, value);
                    }

                    double scaledProbabilitySum = 0.0;
                    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++)
                    {
                        scaledProbabilitySum += Math.Exp(responsibilities[rowIndex, componentIndex] - maximumLogProbability);
                    }
                    if (!IsFinite(scaledProbabilitySum) || scaledProbabilitySum <= 0.0)
                    {
                        throw CreateImpossibleRowException(rowIndex, value);
                    }

                    double rowLogProbability = maximumLogProbability + Math.Log(scaledProbabilitySum);
                    if (!IsFinite(rowLogProbability))
                    {
                        throw CreateImpossibleRowException(rowIndex, value);
                    }

                    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++)
                    {
                        responsibilities[rowIndex, componentIndex] =
                            Math.Exp(responsibilities[rowIndex, componentIndex] - rowLogProbability);
                    }
                    logLikelihood += rowLogProbability;
                }

                return logLikelihood;
            }

            double[] MStep(double[] parameters)
            {
                for (int componentIndex = 0; componentIndex < componentCount; componentIndex++)
                {
                    double weight = 0.0;
                    for (int rowIndex = 0; rowIndex < observationCount; rowIndex++)
                    {
                        if (!IsZeroInflated || sample[rowIndex] > 0.0)
                        {
                            weight += responsibilities[rowIndex, componentIndex];
                        }
                    }
                    mleWeights[componentIndex] = weight;
                }

                double componentWeightSum = mleWeights.Sum();
                double componentWeightTarget = IsZeroInflated ? 1.0 - ZeroWeight : 1.0;
                if (!IsFinite(componentWeightSum) || componentWeightSum <= 0.0)
                {
                    throw new InvalidOperationException("Mixture EM cannot update component weights because no finite positive responsibility mass is available.");
                }

                double scale = componentWeightTarget / componentWeightSum;
                for (int componentIndex = 0; componentIndex < componentCount; componentIndex++)
                {
                    mleWeights[componentIndex] *= scale;
                }

                var solver = new NelderMead(Objective, distributionParameterCount, parameters, lowerParameters, upperParameters);
                solver.Maximize();
                return solver.BestParameterSet.Values;
            }

            double Objective(double[] parameters)
            {
                var distribution = (Mixture)Clone();
                distribution.SetParameters(mleWeights, parameters);
                double logLikelihood = distribution.LogLikelihood(sample);
                return IsFinite(logLikelihood) ? logLikelihood : double.NegativeInfinity;
            }

            InvalidOperationException CreateImpossibleRowException(int rowIndex, double value)
            {
                return new InvalidOperationException(
                    "Mixture EM row " + rowIndex.ToString(CultureInfo.InvariantCulture) +
                    " with value " + value.ToString("R", CultureInfo.InvariantCulture) +
                    " has zero or nonfinite total probability.");
            }

            for (Iterations = 1; Iterations <= MaxIterations; Iterations++)
            {
                newLogLikelihood = EStep(mleParameters);
                if (Math.Abs((oldLogLikelihood - newLogLikelihood) / oldLogLikelihood) < Tolerance) break;
                mleParameters = MStep(mleParameters);
                oldLogLikelihood = newLogLikelihood;
            }

            var result = new List<double>();
            result.AddRange(mleWeights);
            result.AddRange(mleParameters);
            return result.ToArray();
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            if (!_parametersValid) ValidateParameters(GetParameters, true);

            if (IsZeroInflated)
            {
                if (x < 0.0) return 0.0;
                if (x == 0.0) return ZeroWeight;

                double positiveDensity = 0.0;
                for (int i = 0; i < Distributions.Count(); i++)
                {
                    positiveDensity += Weights[i] * PositiveConditionalPDF(i, x);
                }
                return Math.Max(0.0, positiveDensity);
            }

            double density = 0.0;
            for (int i = 0; i < Distributions.Count(); i++) density += Weights[i] * Distributions[i].PDF(x);
            return Math.Max(0.0, density);
        }

        /// <inheritdoc/>
        public override double LogPDF(double x)
        {
            if (!_parametersValid) ValidateParameters(GetParameters, true);

            if (IsZeroInflated)
            {
                if (x < 0.0) return double.NegativeInfinity;
                if (x == 0.0) return Math.Log(ZeroWeight);
            }

            var logDensities = new List<double>();
            for (int i = 0; i < Distributions.Count(); i++)
            {
                double componentLogDensity = IsZeroInflated
                    ? PositiveConditionalLogPDF(i, x)
                    : Distributions[i].LogPDF(x);
                logDensities.Add(Math.Log(Weights[i]) + componentLogDensity);
            }
            return Tools.LogSumExp(logDensities);
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            if (!_parametersValid) ValidateParameters(GetParameters, true);

            if (IsZeroInflated)
            {
                if (x < 0.0) return 0.0;
                if (x == 0.0) return ZeroWeight;

                double hurdleProbability = ZeroWeight;
                for (int i = 0; i < Distributions.Count(); i++)
                {
                    hurdleProbability += Weights[i] * PositiveConditionalCDF(i, x);
                }
                return Clamp(hurdleProbability, 0.0, 1.0);
            }

            double probability = 0.0;
            for (int i = 0; i < Distributions.Count(); i++) probability += Weights[i] * Distributions[i].CDF(x);
            return Clamp(probability, 0.0, 1.0);
        }

        /// <inheritdoc/>
        public override double LogCDF(double x)
        {
            if (!_parametersValid) ValidateParameters(GetParameters, true);
            if (IsZeroInflated) return Math.Log(CDF(x));

            var logProbabilities = new List<double>();
            for (int i = 0; i < Distributions.Count(); i++)
            {
                logProbabilities.Add(Math.Log(Weights[i]) + Distributions[i].LogCDF(x));
            }
            return Tools.LogSumExp(logProbabilities);
        }

        /// <inheritdoc/>
        public override double LogCCDF(double x)
        {
            if (!_parametersValid) ValidateParameters(GetParameters, true);

            if (IsZeroInflated)
            {
                if (x < 0.0) return 0.0;

                double probability = 0.0;
                for (int i = 0; i < Distributions.Count(); i++)
                {
                    probability += Weights[i] * PositiveConditionalCCDF(i, x);
                }
                return Math.Log(Clamp(probability, 0.0, 1.0));
            }

            var logProbabilities = new List<double>();
            for (int i = 0; i < Distributions.Count(); i++)
            {
                logProbabilities.Add(Math.Log(Weights[i]) + Distributions[i].LogCCDF(x));
            }
            return Tools.LogSumExp(logProbabilities);
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            if (probability < 0.0 || probability > 1.0)
                throw new ArgumentOutOfRangeException(nameof(probability), "Probability must be between 0 and 1.");
            if (probability == 0.0) return Minimum;
            if (probability == 1.0) return Maximum;
            if (IsZeroInflated && probability <= ZeroWeight) return 0.0;
            if (!_parametersValid) ValidateParameters(GetParameters, true);

            if (Distributions.Count() == 1 && !IsZeroInflated)
            {
                return Distributions[0].InverseCDF(probability);
            }

            if (_empiricalCDFCreated)
            {
                double empiricalValue = _empiricalCDF.InverseCDF(probability);
                return Clamp(empiricalValue, Minimum, Maximum);
            }

            double componentProbability = IsZeroInflated
                ? (probability - ZeroWeight) / (1.0 - ZeroWeight)
                : probability;
            var componentQuantiles = new List<double>();
            for (int i = 0; i < Distributions.Count(); i++)
            {
                double componentCdfProbability = componentProbability;
                if (IsZeroInflated)
                {
                    componentCdfProbability = PositiveConditionalQuantileProbability(i, componentProbability);
                }
                componentQuantiles.Add(Distributions[i].InverseCDF(componentCdfProbability));
            }

            double lowerBound = componentQuantiles.Min();
            double upperBound = componentQuantiles.Max();
            double value;
            try
            {
                if (lowerBound.AlmostEquals(upperBound)) return Clamp(lowerBound, Minimum, Maximum);
                value = Brent.Solve(y => probability - CDF(y), lowerBound, upperBound, 1E-6, 100, true);
            }
            catch (Exception)
            {
                if (!_empiricalCDFCreated) CreateEmpiricalCDF();
                value = _empiricalCDF.InverseCDF(probability);
            }

            return Clamp(value, Minimum, Maximum);
        }

        /// <summary>
        /// Create empirical distribution for the CDF.
        /// </summary>
        public void CreateEmpiricalCDF()
        {
            // Get min & max
            double minP = 1E-16;
            double maxP = 1 - 1E-16;
            double minX = Distributions.Min(d => d.InverseCDF(minP));
            double maxX = IsZeroInflated
                ? Distributions.Select((distribution, index) =>
                {
                    TryGetPositiveMass(index, out double positiveMass);
                    double probability = distribution.CDF(0.0) + maxP * positiveMass;
                    return distribution.InverseCDF(Math.Min(probability, 1.0 - 1E-15));
                }).Max()
                : Distributions.Max(d => d.InverseCDF(maxP));
            // Get number of bins
            double shift = 0;
            if (minX <= 0) shift = Math.Abs(minX) + 1d;
            double min = minX + shift;
            double max = maxX + shift;
            int order = (int)Math.Floor(Math.Log10(max) - Math.Log10(min));
            int binN = Math.Max(200, 100 * order) - 1;
            // Create bins
            var bins = Stratify.XValues(new StratificationOptions(minX, maxX, binN, false), true);
            var xValues = new List<double>();
            var pValues = new List<double>();
            var x = bins.First().LowerBound;
            var p = CDF(bins.First().LowerBound);
            xValues.Add(x);
            pValues.Add(p);
            for (int i = 1; i < bins.Count; i++)
            {
                x = bins[i].LowerBound;
                p = CDF(x);
                if (x > xValues.Last() && p > pValues.Last())
                {
                    xValues.Add(x);
                    pValues.Add(p);
                }
            }
            x = maxX;
            p = CDF(x);
            if (x > xValues.Last() && p > pValues.Last())
            {
                xValues.Add(x);
                pValues.Add(p);
            }
            _empiricalCDF = new EmpiricalDistribution(xValues, pValues) { XTransform = XTransform, ProbabilityTransform = ProbabilityTransform };
            _empiricalCDFCreated = true;
            _momentsComputed = false;
        }

        /// <inheritdoc/>
        public override double[] GenerateRandomValues(int sampleSize, int seed = -1)
        {
            if (!_parametersValid) ValidateParameters(GetParameters, true);

            var random = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();
            var sample = new double[sampleSize];
            for (int sampleIndex = 0; sampleIndex < sampleSize; sampleIndex++)
            {
                double mixtureProbability = random.NextDouble();
                double componentProbability = random.NextDouble();
                if (IsZeroInflated && mixtureProbability <= ZeroWeight)
                {
                    sample[sampleIndex] = 0.0;
                    continue;
                }

                double cumulativeWeight = IsZeroInflated ? ZeroWeight : 0.0;
                for (int componentIndex = 0; componentIndex < Distributions.Count(); componentIndex++)
                {
                    cumulativeWeight += Weights[componentIndex];
                    if (mixtureProbability <= cumulativeWeight || componentIndex == Distributions.Count() - 1)
                    {
                        double probability = componentProbability;
                        if (IsZeroInflated)
                        {
                            probability = PositiveConditionalQuantileProbability(componentIndex, componentProbability);
                        }
                        sample[sampleIndex] = Distributions[componentIndex].InverseCDF(probability);
                        break;
                    }
                }
            }

            return sample;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            var dists = new UnivariateDistributionBase[Distributions.Count()];
            for (int i = 0; i < Distributions.Count(); i++)
                dists[i] = Distributions[i].Clone();

            return new Mixture(Weights.ToArray(), dists)
            {
                IsZeroInflated = IsZeroInflated,
                ZeroWeight = ZeroWeight,
                XTransform = XTransform,
                ProbabilityTransform = ProbabilityTransform
            };
        }

        /// <inheritdoc/>
        public override XElement ToXElement()
        {
            var result = new XElement("Distribution");
            result.SetAttributeValue(nameof(Type), Type.ToString());
            result.SetAttributeValue(nameof(IsZeroInflated), IsZeroInflated.ToString());
            result.SetAttributeValue(nameof(ZeroWeight), ZeroWeight.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(XTransform), XTransform.ToString());
            result.SetAttributeValue(nameof(ProbabilityTransform), ProbabilityTransform.ToString());
            result.SetAttributeValue(nameof(Distributions), String.Join("|", Distributions.Select(x => x.Type)));
            // Weights
            var weights = Weights;
            var weightStrings = new string[Weights.Length];
            for (int i = 0; i < Weights.Length; i++)
            {
                weightStrings[i] = weights[i].ToString("G17", CultureInfo.InvariantCulture);
            }
            result.SetAttributeValue(nameof(Weights), String.Join("|", weightStrings));
            // Parameters
            var parms = GetParameters;
            var parmStrings = new string[NumberOfParameters];
            for (int i = 0; i < NumberOfParameters; i++)
            {
                parmStrings[i] = parms[i].ToString("G17", CultureInfo.InvariantCulture);
            }
            result.SetAttributeValue("Parameters", String.Join("|", parmStrings));
            return result;
        }

        /// <summary>
        /// Create a mixture distribution from XElement.
        /// </summary>
        /// <param name="xElement">The XElement to deserialize.</param>
        /// <returns>A new mixture distribution.</returns>
        public static Mixture? FromXElement(XElement xElement)
        {
            UnivariateDistributionType type = UnivariateDistributionType.Deterministic;
            var typeAttr = xElement.Attribute(nameof(UnivariateDistributionBase.Type));
            if (typeAttr != null)
            {
                Enum.TryParse(typeAttr.Value, out type);

            }
            if (type == UnivariateDistributionType.Mixture)
            {
                var weights = new List<double>();
                var distributions = new List<UnivariateDistributionBase>();
                var weightsAttr = xElement.Attribute(nameof(Weights));
                if (weightsAttr != null)
                {
                    var w = weightsAttr.Value.Split('|');
                    for (int i = 0; i < w.Length; i++)
                    {
                        double.TryParse(w[i], NumberStyles.Any, CultureInfo.InvariantCulture, out var weight);
                        weights.Add(weight);
                    }
                }
                var distsAttr = xElement.Attribute(nameof(Distributions));
                if (distsAttr != null)
                {
                    var types = distsAttr.Value.Split('|');
                    for (int i = 0; i < types.Length; i++)
                    {
                        Enum.TryParse(types[i], out UnivariateDistributionType distType);
                        distributions.Add(UnivariateDistributionFactory.CreateDistribution(distType));
                    }
                }
                var mixture = new Mixture(weights.ToArray(), distributions.ToArray());

                var zeroInflatedAttr = xElement.Attribute(nameof(IsZeroInflated));
                if (zeroInflatedAttr != null)
                {
                    bool.TryParse(zeroInflatedAttr.Value, out var isZeroInflated);
                    mixture.IsZeroInflated = isZeroInflated;
                }
                var zeroWeightAttr = xElement.Attribute(nameof(ZeroWeight));
                if (zeroWeightAttr != null)
                {
                    double.TryParse(zeroWeightAttr.Value, NumberStyles.Any, CultureInfo.InvariantCulture, out var zeroWeight);
                    mixture.ZeroWeight = zeroWeight;
                }
                var xTransformAttr = xElement.Attribute(nameof(XTransform));
                if (xTransformAttr != null)
                {
                    Enum.TryParse(xTransformAttr.Value, out Transform xTransform);
                    mixture.XTransform = xTransform;
                }
                var probTransformAttr = xElement.Attribute(nameof(ProbabilityTransform));
                if (probTransformAttr != null)
                {
                    Enum.TryParse(probTransformAttr.Value, out Transform probabilityTransform);
                    mixture.ProbabilityTransform = probabilityTransform;
                }
                var paramsAttr = xElement.Attribute("Parameters");
                if (paramsAttr != null)
                {
                    var vals = paramsAttr.Value.Split('|');
                    var parameters = new List<double>();
                    for (int i = 0; i < vals.Length; i++)
                    {
                        double.TryParse(vals[i], NumberStyles.Any, CultureInfo.InvariantCulture, out var parm);
                        parameters.Add(parm);
                    }
                    mixture.SetParameters(parameters);
                }

                return mixture;
            }
            else
            {
                return null;
            }
        }

    }
}
