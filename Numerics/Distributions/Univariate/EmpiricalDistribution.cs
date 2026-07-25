using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Globalization;
using System.Linq;
using System.Xml.Linq;
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
    /// a probability. Points on the cumulative curve must be entered with nondecreasing values and increasing
    /// probabilities. Even though the (X,p) points define the distribution, any value between the minimum
    /// and maximum can be returned.
    /// </para>
    /// <para>
    /// <b> References:</b>
    /// <list type="bullet">
    /// <item><description>
    /// The distribution behaves similarly to the "RiskCumul" function in the Palisade's @Risk software.
    /// <see href="http://kb.palisade.com/index.php?pg=kb.page&amp;id=51"/>
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
        /// <param name="XOrder">Ascending sort order of X values.</param>
        /// <param name="probabilityOrder">Sort order of probability values.</param>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when the distribution is used and <paramref name="XOrder"/> is not ascending.
        /// </exception>
        public EmpiricalDistribution(IList<double> XValues, IList<double> PValues, SortOrder XOrder, SortOrder probabilityOrder)
        {
            SetParameters(XValues, PValues, XOrder, probabilityOrder);

        }

        /// <summary>
        /// Constructs a Univariate Empirical CDF from ordered paired data.
        /// </summary>
        /// <param name="orderedPairedData">The ordered paired data.</param>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when the X values are not configured in ascending order.
        /// </exception>
        public EmpiricalDistribution(OrderedPairedData orderedPairedData)
        {
            if (orderedPairedData.OrderX != SortOrder.Ascending) throw new ArgumentOutOfRangeException(nameof(orderedPairedData), "Empirical x-values must be ordered ascending.");
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
                    opd = new OrderedPairedData(_xValues, _pValues, orderedPairedData.StrictX, SortOrder.Ascending, true, SortOrder.Ascending);
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
            _parametersValid = ValidateData(opd, _pValues, false) is null;
            _momentsComputed = false;
        }

        /// <summary>
        /// Constructs a Univariate Empirical CDF from sample data.
        /// </summary>
        /// <param name="sample">The sample values used as empirical ordinates; duplicate values are permitted.</param>
        /// <param name="plottingPostionType">The plotting position formula type. Default = Weibull.</param>
        public EmpiricalDistribution(IList<double> sample, PlottingPositions.PlottingPostionType plottingPostionType = PlottingPositions.PlottingPostionType.Weibull)
        {
            _xValues = sample.ToArray();
            Array.Sort(_xValues);
            _pValues = PlottingPositions.Function(_xValues.Count(), plottingPostionType)!;
            opd = new OrderedPairedData(_xValues, _pValues, false, SortOrder.Ascending, true, SortOrder.Ascending);
            _parametersValid = ValidateData(opd, _pValues, false) is null;
            _momentsComputed = false;
        }

        private double[] _xValues = null!;
        private double[] _pValues = null!;
        private OrderedPairedData opd = null!;
        private bool _momentsComputed = false;
        private double u1, u2, u3, u4;

        /// <summary>
        /// Returns the array of X values. Points On the cumulative curve are specified
        /// with nondecreasing values and increasing probabilities.
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
            get { return []; }
        }


        /// <summary>
        /// Compute central moments of the distribution.
        /// </summary>
        private void ComputeMoments()
        {
            var mom = CentralMoments(300);
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
        /// <exception cref="ArgumentException">The value and probability collections have different lengths.</exception>
        public void SetParameters(IList<double> xValues, IList<double> pValues)
        {
            if (xValues.Count != pValues.Count)
            {
                throw new ArgumentException("The value and probability arrays must have the same length.", nameof(pValues));
            }
            _xValues = xValues.ToArray();
            _pValues = pValues.ToArray();
            opd = new OrderedPairedData(xValues, pValues, false, SortOrder.Ascending, true, SortOrder.Ascending);
            _parametersValid = ValidateData(opd, pValues, false) is null;
            _momentsComputed = false;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="xValues">Array of X values.</param>
        /// <param name="pValues">Array of probability values. Range 0 ≤ p ≤ 1.</param>
        /// <param name="XOrder">Ascending sort order of X values.</param>
        /// <param name="probabilityOrder">Sort order of probability values.</param>
        /// <exception cref="ArgumentException">The value and probability collections have different lengths.</exception>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when the distribution is used and <paramref name="XOrder"/> is not ascending.
        /// </exception>
        public void SetParameters(IList<double> xValues, IList<double> pValues, SortOrder XOrder, SortOrder probabilityOrder)
        {
            if (xValues.Count != pValues.Count)
            {
                throw new ArgumentException("The value and probability arrays must have the same length.", nameof(pValues));
            }
            _xValues = xValues.ToArray();
            _pValues = pValues.ToArray();
            opd = new OrderedPairedData(xValues, pValues, false, XOrder, true, probabilityOrder);
            _parametersValid = ValidateData(opd, pValues, false) is null;
            _momentsComputed = false;
        }

        /// <summary>
        /// Validates empirical ordinates and their cumulative probabilities.
        /// </summary>
        /// <param name="data">The ordered paired data to validate.</param>
        /// <param name="probabilities">The cumulative probabilities.</param>
        /// <param name="throwException">Whether to throw the validation exception.</param>
        /// <returns><see langword="null"/> when the data are valid; otherwise, the validation exception.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="throwException"/> is true and the data are invalid.</exception>
        private static ArgumentOutOfRangeException? ValidateData(OrderedPairedData data, IList<double> probabilities, bool throwException)
        {
            ArgumentOutOfRangeException? exception = null;
            if (data.Count < 2)
            {
                exception = new ArgumentOutOfRangeException(nameof(data), "At least two empirical ordinates are required.");
            }
            else if (probabilities.Count != data.Count)
            {
                exception = new ArgumentOutOfRangeException(nameof(probabilities), "The empirical ordinate and probability collections must have the same length.");
            }
            else if (data.OrderX != SortOrder.Ascending)
            {
                exception = new ArgumentOutOfRangeException(nameof(data), "Empirical x-values must be ordered ascending.");
            }
            else if (!data.IsValid)
            {
                List<string> errors = data.GetErrors().Distinct().ToList();
                string message = errors.Count > 0
                    ? string.Join(" ", errors)
                    : "The empirical ordinates do not satisfy their configured ordering.";
                exception = new ArgumentOutOfRangeException(nameof(data), message);
            }
            else if (probabilities.Any(probability =>
                double.IsNaN(probability) ||
                double.IsInfinity(probability) ||
                probability < 0d ||
                probability > 1d))
            {
                exception = new ArgumentOutOfRangeException(nameof(probabilities), "Empirical probabilities must be finite and between zero and one.");
            }
            if (throwException && exception is not null) throw exception;
            return exception;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            throw new NotImplementedException();
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
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
            if (_parametersValid == false) ValidateData(opd, _pValues, true);
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
            if (_parametersValid == false) ValidateData(opd, _pValues, true);
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

        /// <summary>
        /// Convolves two empirical distributions with a log-spaced output grid — the opt-in
        /// alternative for order-of-magnitude (heavy-tail) supports that the linear output grid
        /// under-resolves. The FFT pipeline is exactly the linear-grid
        /// <see cref="Convolve(EmpiricalDistribution, EmpiricalDistribution, int)"/> (unchanged,
        /// byte-identical when this overload is not used); only the OUTPUT resampling differs:
        /// the returned distribution's CDF is re-read on a log-spaced ladder between the summed
        /// supports, keeping tail resolution across decades.
        /// </summary>
        /// <param name="dist1">The first empirical distribution.</param>
        /// <param name="dist2">The second empirical distribution.</param>
        /// <param name="numberOfPoints">The number of output points. Default = 1024.</param>
        /// <param name="logSpacedOutput">True for the log-spaced output ladder; false delegates to the linear-grid method unchanged.</param>
        /// <returns>The convolved empirical distribution.</returns>
        /// <exception cref="ArgumentException">Thrown when a log-spaced output is requested over a non-positive summed support.</exception>
        public static EmpiricalDistribution Convolve(EmpiricalDistribution dist1, EmpiricalDistribution dist2, int numberOfPoints, bool logSpacedOutput)
        {
            var linear = Convolve(dist1, dist2, numberOfPoints);
            if (!logSpacedOutput) return linear;

            double minimum = linear.Minimum;
            double maximum = linear.Maximum;
            if (minimum <= 0d || maximum <= 0d)
                throw new ArgumentException("A log-spaced output grid requires a strictly positive summed support.", nameof(logSpacedOutput));

            double logMin = Math.Log10(minimum);
            double logMax = Math.Log10(maximum);
            var xValues = new double[numberOfPoints];
            var pValues = new double[numberOfPoints];
            double previous = double.NegativeInfinity;
            int count = 0;
            for (int i = 0; i < numberOfPoints; i++)
            {
                double x = Math.Pow(10d, logMin + (logMax - logMin) * i / (numberOfPoints - 1d));
                double p = linear.CDF(x);
                // The CDF ladder must stay strictly increasing for the empirical constructor.
                if (p > previous)
                {
                    xValues[count] = x;
                    pValues[count] = p;
                    previous = p;
                    count++;
                }
            }
            var trimmedX = new double[count];
            var trimmedP = new double[count];
            Array.Copy(xValues, trimmedX, count);
            Array.Copy(pValues, trimmedP, count);
            return new EmpiricalDistribution(trimmedX, trimmedP) { XTransform = dist1.XTransform, ProbabilityTransform = dist1.ProbabilityTransform };
        }

        /// <summary>
        /// Convolves two discrete (atom-bearing) distributions EXACTLY on a shared uniform
        /// lattice — the atom-aware form the continuous <c>Convolve</c> cannot represent: it
        /// samples continuous PDFs, so a point mass (a CDF jump, e.g. the zero-inflation atom of
        /// a defective risk curve) has no representation there. Each input's atoms deposit onto
        /// the lattice with a moment-preserving two-node split (input means are preserved
        /// exactly), the mass vectors convolve by FFT, and the result is the discrete mass
        /// ladder on the summed lattice.
        /// </summary>
        /// <param name="values1">The first distribution's atom values.</param>
        /// <param name="masses1">The first distribution's atom masses (non-negative; typically summing to one).</param>
        /// <param name="values2">The second distribution's atom values.</param>
        /// <param name="masses2">The second distribution's atom masses.</param>
        /// <param name="latticePoints">The per-input lattice resolution (rounded up to a power of two). Default = 4096.</param>
        /// <param name="values">The convolved lattice values (uniform, ascending).</param>
        /// <param name="masses">The convolved lattice masses (non-negative; FFT ringing is floored and the total mass renormalized to the product of the input totals).</param>
        /// <exception cref="ArgumentNullException">Thrown when any input list is null.</exception>
        /// <exception cref="ArgumentException">Thrown when a values/masses pair is empty or mismatched in length, or a mass is negative or non-finite.</exception>
        public static void ConvolveDiscrete(IList<double> values1, IList<double> masses1, IList<double> values2, IList<double> masses2,
            int latticePoints, out double[] values, out double[] masses)
        {
            ValidateAtoms(values1, masses1, nameof(values1));
            ValidateAtoms(values2, masses2, nameof(values2));
            if (latticePoints < 8) latticePoints = 8;
            int n = Tools.NextPowerOfTwo(latticePoints);

            // The shared lattice: each input is binned over its own span, but on the SAME step
            // so the convolution lattice is uniform. The step spans the summed support.
            double min1 = Min(values1), max1 = Max(values1);
            double min2 = Min(values2), max2 = Max(values2);
            double span = Math.Max((max1 - min1) + (max2 - min2), Tools.DoubleMachineEpsilon);
            double step = span / (n - 1);

            var lattice1 = DepositAtoms(values1, masses1, min1, step, n);
            var lattice2 = DepositAtoms(values2, masses2, min2, step, n);

            // FFT convolution of the mass vectors (zero-padded complex arrays).
            int fftSize = Tools.NextPowerOfTwo(2 * n);
            var fft1 = new double[2 * fftSize];
            var fft2 = new double[2 * fftSize];
            for (int i = 0; i < n; i++)
            {
                fft1[2 * i] = lattice1[i];
                fft2[2 * i] = lattice2[i];
            }
            Mathematics.Fourier.FFT(fft1);
            Mathematics.Fourier.FFT(fft2);
            var product = new double[2 * fftSize];
            for (int i = 0; i < fftSize; i++)
            {
                double re1 = fft1[2 * i], im1 = fft1[2 * i + 1];
                double re2 = fft2[2 * i], im2 = fft2[2 * i + 1];
                product[2 * i] = re1 * re2 - im1 * im2;
                product[2 * i + 1] = re1 * im2 + im1 * re2;
            }
            Mathematics.Fourier.FFT(product, inverse: true);

            // The result ladder: 2n − 1 meaningful nodes from min1 + min2, floored against FFT
            // ringing and renormalized to the exact product of the input mass totals.
            int resultCount = 2 * n - 1;
            values = new double[resultCount];
            masses = new double[resultCount];
            double total = 0d;
            for (int i = 0; i < resultCount; i++)
            {
                values[i] = (min1 + min2) + i * step;
                double mass = product[2 * i] / fftSize;
                masses[i] = mass > 0d ? mass : 0d;
                total += masses[i];
            }
            double expected = Sum(masses1) * Sum(masses2);
            if (total > 0d && expected > 0d)
            {
                double scale = expected / total;
                for (int i = 0; i < resultCount; i++) masses[i] *= scale;
            }
        }

        /// <summary>
        /// Validates one atom list pair.
        /// </summary>
        /// <param name="values">The atom values.</param>
        /// <param name="masses">The atom masses.</param>
        /// <param name="parameterName">The reported parameter name.</param>
        private static void ValidateAtoms(IList<double> values, IList<double> masses, string parameterName)
        {
            if (values == null || masses == null) throw new ArgumentNullException(parameterName);
            if (values.Count == 0 || values.Count != masses.Count)
                throw new ArgumentException("Atom values and masses must be non-empty and equal in length.", parameterName);
            for (int i = 0; i < masses.Count; i++)
            {
                if (!Tools.IsFinite(values[i]) || !Tools.IsFinite(masses[i]) || masses[i] < 0d)
                    throw new ArgumentException("Atom values must be finite and masses non-negative.", parameterName);
            }
        }

        /// <summary>
        /// Deposits atoms onto a uniform lattice with the moment-preserving two-node split: an
        /// atom between nodes splits its mass so the lattice mean reproduces the atom mean
        /// exactly.
        /// </summary>
        /// <param name="values">The atom values.</param>
        /// <param name="masses">The atom masses.</param>
        /// <param name="origin">The lattice origin.</param>
        /// <param name="step">The lattice step.</param>
        /// <param name="n">The lattice node count.</param>
        /// <returns>The lattice mass vector.</returns>
        private static double[] DepositAtoms(IList<double> values, IList<double> masses, double origin, double step, int n)
        {
            var lattice = new double[n];
            for (int i = 0; i < values.Count; i++)
            {
                double position = (values[i] - origin) / step;
                int lower = (int)Math.Floor(position);
                if (lower < 0) lower = 0;
                if (lower > n - 1) lower = n - 1;
                int upper = lower + 1;
                if (upper > n - 1)
                {
                    lattice[n - 1] += masses[i];
                    continue;
                }
                double fraction = position - lower;
                if (fraction < 0d) fraction = 0d;
                if (fraction > 1d) fraction = 1d;
                lattice[lower] += masses[i] * (1d - fraction);
                lattice[upper] += masses[i] * fraction;
            }
            return lattice;
        }

        /// <summary>Returns the smallest value of a list.</summary>
        /// <param name="values">The list.</param>
        private static double Min(IList<double> values)
        {
            double minimum = double.MaxValue;
            for (int i = 0; i < values.Count; i++) if (values[i] < minimum) minimum = values[i];
            return minimum;
        }

        /// <summary>Returns the largest value of a list.</summary>
        /// <param name="values">The list.</param>
        private static double Max(IList<double> values)
        {
            double maximum = double.MinValue;
            for (int i = 0; i < values.Count; i++) if (values[i] > maximum) maximum = values[i];
            return maximum;
        }

        /// <summary>Returns the sum of a list.</summary>
        /// <param name="values">The list.</param>
        private static double Sum(IList<double> values)
        {
            double sum = 0d;
            for (int i = 0; i < values.Count; i++) sum += values[i];
            return sum;
        }

        /// <summary>
        /// Serializes the empirical distribution to an XElement, overriding the base scalar-only
        /// form: the X and probability tables (round-trip-exact "G17"), the probability sort
        /// order, and both interpolation transforms. The base implementation writes only scalar
        /// parameters, which lost the tables entirely.
        /// </summary>
        /// <returns>An XElement representation of the empirical distribution.</returns>
        public override XElement ToXElement()
        {
            var result = new XElement("Distribution");
            result.SetAttributeValue(nameof(Type), Type.ToString());
            result.SetAttributeValue(nameof(XTransform), XTransform.ToString());
            result.SetAttributeValue(nameof(ProbabilityTransform), ProbabilityTransform.ToString());

            var xValues = new string[XValues.Count];
            var pValues = new string[ProbabilityValues.Count];
            for (int i = 0; i < XValues.Count; i++)
                xValues[i] = XValues[i].ToString("G17", CultureInfo.InvariantCulture);
            for (int i = 0; i < ProbabilityValues.Count; i++)
                pValues[i] = ProbabilityValues[i].ToString("G17", CultureInfo.InvariantCulture);
            result.SetAttributeValue(nameof(XValues), string.Join("|", xValues));
            result.SetAttributeValue(nameof(ProbabilityValues), string.Join("|", pValues));

            // The stored probability ladder may run ascending (non-exceedance) or descending
            // (exceedance); record the order so deserialization restores the same convention.
            var order = SortOrder.Ascending;
            if (ProbabilityValues.Count > 1 && ProbabilityValues[0] > ProbabilityValues[ProbabilityValues.Count - 1])
                order = SortOrder.Descending;
            result.SetAttributeValue("ProbabilityOrder", order.ToString());
            return result;
        }

        /// <summary>
        /// Deserializes an empirical distribution from an XElement produced by
        /// <see cref="ToXElement"/>.
        /// </summary>
        /// <param name="xElement">The XElement to deserialize.</param>
        /// <returns>A new <see cref="EmpiricalDistribution"/>.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="xElement"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the element carries no parseable X/probability tables.</exception>
        public static EmpiricalDistribution FromXElement(XElement xElement)
        {
            if (xElement == null) throw new ArgumentNullException(nameof(xElement));
            string? xText = xElement.Attribute(nameof(XValues))?.Value;
            string? pText = xElement.Attribute(nameof(ProbabilityValues))?.Value;
            if (string.IsNullOrEmpty(xText) || string.IsNullOrEmpty(pText))
                throw new ArgumentException("The serialized empirical distribution is missing its X or probability table.", nameof(xElement));

            string[] xTokens = xText!.Split('|');
            string[] pTokens = pText!.Split('|');
            if (xTokens.Length != pTokens.Length)
                throw new ArgumentException("The serialized empirical distribution's X and probability tables differ in length.", nameof(xElement));

            var xValues = new double[xTokens.Length];
            var pValues = new double[pTokens.Length];
            for (int i = 0; i < xTokens.Length; i++)
            {
                if (!double.TryParse(xTokens[i], NumberStyles.Any, CultureInfo.InvariantCulture, out xValues[i])
                    || !double.TryParse(pTokens[i], NumberStyles.Any, CultureInfo.InvariantCulture, out pValues[i]))
                {
                    throw new ArgumentException("The serialized empirical distribution carries an unparseable table value.", nameof(xElement));
                }
            }

            var order = SortOrder.Ascending;
            var orderAttr = xElement.Attribute("ProbabilityOrder");
            if (orderAttr != null) Enum.TryParse(orderAttr.Value, out order);

            var distribution = new EmpiricalDistribution(xValues, pValues, SortOrder.Ascending, order);
            if (Enum.TryParse(xElement.Attribute(nameof(XTransform))?.Value, out Transform xTransform))
                distribution.XTransform = xTransform;
            if (Enum.TryParse(xElement.Attribute(nameof(ProbabilityTransform))?.Value, out Transform probabilityTransform))
                distribution.ProbabilityTransform = probabilityTransform;
            return distribution;
        }


        /// <summary>
        /// Convolves two empirical distributions using FFT.
        /// </summary>
        /// <param name="dist1">The first empirical distribution.</param>
        /// <param name="dist2">The second empirical distribution.</param>
        /// <param name="numberOfPoints">Number of points in the resulting distribution. Default = 1024.</param>
        /// <returns>A new empirical distribution representing the convolution.</returns>
        public static EmpiricalDistribution Convolve(EmpiricalDistribution dist1, EmpiricalDistribution dist2, int numberOfPoints = 1024)
        {
            if (numberOfPoints < 2)
                throw new ArgumentException("Number of points must be at least 2.", nameof(numberOfPoints));

            int requestedPoints = numberOfPoints;

            // Use a large FFT size for accuracy - at least 8x the requested points
            int fftPoints = Tools.NextPowerOfTwo(Math.Max(numberOfPoints * 8, 2048));

            // Result range
            double minResult = dist1.Minimum + dist2.Minimum;
            double maxResult = dist1.Maximum + dist2.Maximum;
            double rangeResult = maxResult - minResult;

            // Create uniform grid over result range
            double dx = rangeResult / fftPoints;

            // Sample both PDFs on a uniform grid
            // For FFT convolution, we need both PDFs on the same spacing
            var pdf1 = new double[fftPoints];
            var pdf2 = new double[fftPoints];

            // Sample each PDF on its own domain, then pad with zeros
            for (int i = 0; i < fftPoints; i++)
            {
                double x1 = dist1.Minimum + i * dx;
                double x2 = dist2.Minimum + i * dx;

                if (x1 >= dist1.Minimum && x1 <= dist1.Maximum)
                {
                    pdf1[i] = dist1.PDF(x1);
                }

                if (x2 >= dist2.Minimum && x2 <= dist2.Maximum)
                {
                    pdf2[i] = dist2.PDF(x2);
                }
            }

            // Verify and normalize PDFs
            double integral1 = 0, integral2 = 0;
            for (int i = 0; i < fftPoints - 1; i++)
            {
                integral1 += 0.5 * (pdf1[i] + pdf1[i + 1]) * dx;
                integral2 += 0.5 * (pdf2[i] + pdf2[i + 1]) * dx;
            }

            // Normalize if needed
            if (integral1 > 0.01)
            {
                for (int i = 0; i < fftPoints; i++)
                {
                    pdf1[i] /= integral1;
                }
            }
            if (integral2 > 0.01)
            {
                for (int i = 0; i < fftPoints; i++)
                {
                    pdf2[i] /= integral2;
                }
            }

            // Prepare for FFT with zero padding to avoid circular convolution
            int fftSize = Tools.NextPowerOfTwo(fftPoints * 2);
            var fft1 = new double[2 * fftSize];
            var fft2 = new double[2 * fftSize];

            for (int i = 0; i < fftPoints; i++)
            {
                fft1[2 * i] = pdf1[i] * dx; // Scale by dx for discrete convolution
                fft2[2 * i] = pdf2[i] * dx;
            }

            // Perform FFT
            Fourier.FFT(fft1);
            Fourier.FFT(fft2);

            // Multiply in frequency domain
            var result = new double[2 * fftSize];
            for (int i = 0; i < fftSize; i++)
            {
                int idx = 2 * i;
                double real1 = fft1[idx];
                double imag1 = fft1[idx + 1];
                double real2 = fft2[idx];
                double imag2 = fft2[idx + 1];

                result[idx] = real1 * real2 - imag1 * imag2;
                result[idx + 1] = real1 * imag2 + imag1 * real2;
            }

            // Inverse FFT
            Fourier.FFT(result, inverse: true);

            // Extract and normalize convolved PDF
            var convolvedPDF = new double[fftSize];
            double maxPDF = 0;
            for (int i = 0; i < fftSize; i++)
            {
                convolvedPDF[i] = result[2 * i] / fftSize;
                if (convolvedPDF[i] < 0) convolvedPDF[i] = 0;
                if (convolvedPDF[i] > maxPDF) maxPDF = convolvedPDF[i];
            }

            // Build X-axis for convolved PDF - starts at the minimum of the result range
            var xConv = new double[fftSize];
            double dxFine = dx; // Use same spacing as input
            for (int i = 0; i < fftSize; i++)
            {
                xConv[i] = dist1.Minimum + dist2.Minimum + i * dxFine;
            }

            // Integrate to get CDF
            var cdfFine = new double[fftSize];
            cdfFine[0] = 0;

            for (int i = 1; i < fftSize; i++)
            {
                cdfFine[i] = cdfFine[i - 1] + 0.5 * (convolvedPDF[i - 1] + convolvedPDF[i]) * dxFine;
            }

            // Normalize CDF
            double cdfMax = cdfFine[fftSize - 1];
            if (cdfMax > 1E-10)
            {
                for (int i = 0; i < fftSize; i++)
                {
                    cdfFine[i] /= cdfMax;
                }
            }

            // Resample to requested resolution using interpolation
            var finalXValues = new double[requestedPoints];
            var finalCdfValues = new double[requestedPoints];

            // Only use points within the valid range for interpolation
            var validX = new List<double>();
            var validCDF = new List<double>();

            for (int i = 0; i < fftSize; i++)
            {
                if (xConv[i] >= minResult && xConv[i] <= maxResult)
                {
                    validX.Add(xConv[i]);
                    validCDF.Add(cdfFine[i]);
                }
            }

            var tempOpd = new OrderedPairedData(
                validX.ToArray(),
                validCDF.ToArray(),
                true, SortOrder.Ascending,
                true, SortOrder.Ascending);

            double outputDx = (maxResult - minResult) / (requestedPoints - 1);
            for (int i = 0; i < requestedPoints; i++)
            {
                finalXValues[i] = minResult + i * outputDx;
                finalCdfValues[i] = tempOpd.GetYFromX(finalXValues[i], Transform.None, Transform.None);

                // Ensure valid CDF values
                if (finalCdfValues[i] < 0) finalCdfValues[i] = 0;
                if (finalCdfValues[i] > 1) finalCdfValues[i] = 1;
            }

            // Ensure monotonicity
            for (int i = 1; i < requestedPoints; i++)
            {
                if (finalCdfValues[i] < finalCdfValues[i - 1])
                {
                    finalCdfValues[i] = finalCdfValues[i - 1];
                }
            }

            return new EmpiricalDistribution(finalXValues, finalCdfValues, SortOrder.Ascending, SortOrder.Ascending);
        }

        /// <summary>
        /// Convolves a list of empirical distributions using FFT.
        /// </summary>
        /// <param name="distributions">List of empirical distributions to convolve.</param>
        /// <param name="numberOfPoints">Number of points in the resulting distribution. Default = 1024.</param>
        /// <returns>A new empirical distribution representing the convolution of all distributions.</returns>
        public static EmpiricalDistribution Convolve(IList<EmpiricalDistribution> distributions, int numberOfPoints = 1024)
        {
            if (distributions == null || distributions.Count == 0)
                throw new ArgumentException("Distribution list cannot be null or empty.", nameof(distributions));

            if (distributions.Count == 1)
                return (EmpiricalDistribution)distributions[0].Clone();

            if (numberOfPoints < 2)
                throw new ArgumentException("Number of points must be at least 2.", nameof(numberOfPoints));

            // Start with the first distribution
            var result = distributions[0];

            // Convolve sequentially with each subsequent distribution
            for (int i = 1; i < distributions.Count; i++)
            {
                result = Convolve(result, distributions[i], numberOfPoints);
            }

            return result;
        }

    }
}