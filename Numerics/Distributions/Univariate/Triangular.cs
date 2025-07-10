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

using Numerics.Mathematics.Optimization;
using System;
using System.Collections.Generic;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Triangular probability distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> References:</b>
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Triangular_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class Triangular : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IBootstrappable
    {
       
        /// <summary>
        /// Constructs a Triangular distribution with min = 0.0, max = 1.0, and mode = 0.5.
        /// </summary>
        public Triangular()
        {
            SetParameters(0.0d, 0.5d, 1.0d);
        }

        /// <summary>
        /// Constructs a Triangular distribution with specified min, max, and mode.
        /// </summary>
        /// <param name="min">The minimum possible value of the distribution.</param>
        /// <param name="mode">The mode, or most likely, value of the distribution.</param>
        /// <param name="max">The maximum possible value of the distribution.</param>
        /// <remarks>
        /// In probability theory and statistics, the triangular distribution is a continuous probability distribution
        /// with lower limit a, upper limit b and mode c.
        /// </remarks>
        public Triangular(double min, double mode, double max)
        {
            // Set parameters
            SetParameters(min, mode, max);
        }

        private double _min;
        private double _max;
        private double _mode;

        /// <summary>
        /// Get and set the min of the distribution.
        /// </summary>
        public double Min
        {
            get { return _min; }
            set
            {
                _parametersValid = ValidateParameters(value, MostLikely, Max, false) is null;
                _min = value;
            }
        }

        /// <summary>
        /// Get and set the mode of the distribution.
        /// </summary>
        public double MostLikely
        {
            get { return _mode; }
            set
            {
                _parametersValid = ValidateParameters(Min, value, Max, false) is null;
                _mode = value;
            }
        }

        /// <summary>
        /// Get and set the max of the distribution.
        /// </summary>
        public double Max
        {
            get { return _max; }
            set
            {
                _parametersValid = ValidateParameters(Min, MostLikely, value, false) is null;
                _max = value;
            }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 3; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.Triangular; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Triangular"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "TRI"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[3, 2];
                parmString[0, 0] = "Min (a)";
                parmString[1, 0] = "Most Likely (c)";
                parmString[2, 0] = "Max (b)";
                parmString[0, 1] = Min.ToString();
                parmString[1, 1] = MostLikely.ToString();
                parmString[2, 1] = Max.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["a", "c", "b"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Min), nameof(MostLikely), nameof(Max)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Min, MostLikely, Max]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get { return (_min + _max + _mode) / 3.0d; }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get
            {
                double m = 0;
                if (_mode >= (_min + _max) / 2.0d)
                {
                    m = _min + Math.Sqrt((_max - _min) * (_mode - _min) / 2.0d);
                }
                else if (_mode <= (_min + _max) / 2.0d)
                {
                    m = _max - Math.Sqrt((_max - _min) * (_max - _mode) / 2.0d);
                }
                return m;
            }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get { return _mode; }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get { return Math.Sqrt((_min * _min + _max * _max + _mode * _mode - _min * _max - _min * _mode - _max * _mode) / 18.0d); }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {
                double q = Tools.Sqrt2 * (_min + _max - 2d * _mode) * (2d * _min - _max - _mode) * (_min - 2d * _max + _mode);
                double d = 5d * Math.Pow(_min * _min + _max * _max + _mode * _mode - _min * _max - _min * _mode - _max * _mode, 3.0d / 2.0d);
                return q / d;
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get { return 12d / 5d; }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return _min; }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return _max; }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [double.NegativeInfinity, Min, Mode]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [Mode, Max, double.PositiveInfinity]; }
        }

        /// <inheritdoc/>
        public void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod)
        {
            if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
            {
                var min = Tools.Min(sample);
                var max = Tools.Max(sample);
                var mean = Tools.Mean(sample);
                var mode = mean * 3 - max - min;
                mode = Math.Max(min, mode);
                mode = Math.Min(max, mode);
                SetParameters(min, mode, max);
            }
            else if (estimationMethod == ParameterEstimationMethod.MaximumLikelihood)
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
            var newDistribution = (Triangular)Clone();
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="min">The minimum possible value of the distribution.</param>
        /// <param name="mode">The mode, or most likely, value of the distribution.</param>
        /// <param name="max">The maximum possible value of the distribution.</param>
        public void SetParameters(double min, double mode, double max)
        {
            // validate parameters
            _parametersValid = ValidateParameters(min, mode, max, false) is null;
            // Set parameters
            _min = min;
            _mode = mode;
            _max = max;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1], parameters[2]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="min">The minimum possible value of the distribution.</param>
        /// <param name="mode">The mode, or most likely, value of the distribution.</param>
        /// <param name="max">The maximum possible value of the distribution.</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        private ArgumentOutOfRangeException ValidateParameters(double min, double mode, double max, bool throwException)
        {
            if (double.IsNaN(min) || double.IsInfinity(min) ||
                double.IsNaN(max) || double.IsInfinity(max) || min > max)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Min), "The min cannot be greater than the max.");
                return new ArgumentOutOfRangeException(nameof(Min), "The min cannot be greater than the max.");
            }
            else if (double.IsNaN(mode) || double.IsInfinity(mode) || mode < min || mode > max)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(MostLikely), "The mode (most likely) must be between the min and max.");
                return new ArgumentOutOfRangeException(nameof(MostLikely), "The mode (most likely) must be between the min and max.");
            }
            return null;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            return ValidateParameters(parameters[0], parameters[1], parameters[2], throwException);
        }

        /// <inheritdoc/>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            var initialVals = new double[NumberOfParameters];
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];
            // Estimate initial values using the method of moments (a.k.a product moments).
            var min = Tools.Min(sample);
            var max = Tools.Max(sample);
            var mean = Tools.Mean(sample);
            var mode = mean * 3 - max - min;
            initialVals[0] = min - Tools.DoubleMachineEpsilon;
            initialVals[1] = mode;
            initialVals[2] = max + Tools.DoubleMachineEpsilon;
            // Get bounds of min
            lowerVals[0] = -Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(min)) + 1d));
            upperVals[0] = min;
            // Get bounds of mode
            lowerVals[1] = min;
            upperVals[1] = max;
            // Get bounds of max
            lowerVals[2] = max;
            upperVals[2] = Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(max)) + 1d));

            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
        }

        /// <inheritdoc/>
        public double[] MLE(IList<double> sample)
        {
            // Set constraints
            var tuple = GetParameterConstraints(sample);
            var Initials = tuple.Item1;
            var Lowers = tuple.Item2;
            var Uppers = tuple.Item3;

            // Solve using Nelder-Mead (Downhill Simplex)
            double logLH(double[] x)
            {
                var N = new Triangular();
                N.SetParameters(x);
                return N.LogLikelihood(sample);
            }
            var solver = new NelderMead(logLH, NumberOfParameters, Initials, Lowers, Uppers);
            solver.ReportFailure = true;
            solver.Maximize();
            return solver.BestParameterSet.Values;
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Min, MostLikely, Max, true);
            // 
            if (Min.AlmostEquals(Max) && Min.AlmostEquals(MostLikely)) return 0.0d;
            if (x < Minimum || x > Maximum) return 0.0d;
            if (x >= Min && x < MostLikely)
            {
                return 2.0d * (x - Min) / ((Max - Min) * (MostLikely - Min));
            }
            else if (x > MostLikely && x <= Max)
            {
                return 2.0d * (Max - x) / ((Max - Min) * (Max - MostLikely));
            }
            else if (x == MostLikely)
            {
                return 2.0d / (Max - Min);
            }
            return double.NaN;
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Min, MostLikely, Max, true);
            if (Min.AlmostEquals(Max) && Min.AlmostEquals(MostLikely)) return 1d;         
            if (x <= Minimum) return 0d;
            if (x >= Maximum) return 1d;
            if (x > Min && x <= MostLikely)
            {
                return (x - Min) * (x - Min) / ((Max - Min) * (MostLikely - Min));
            }
            else if (x > MostLikely && x < Max)
            {
                return 1d - (Max - x) * (Max - x) / ((Max - Min) * (Max - MostLikely));
            }
            return double.NaN;
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (probability < 0.0d || probability > 1.0d)
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            if (Min.AlmostEquals(Max) && Min.AlmostEquals(MostLikely)) return Min;
            if (probability == 0.0d) return Minimum;
            if (probability == 1.0d) return Maximum;
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Min, MostLikely, Max, true);
            if (probability < (MostLikely - Min) / (Max - Min))
            {
                return Min + Math.Sqrt(probability * (Max - Min) * (MostLikely - Min));
            }
            else if (probability >= (MostLikely - Min) / (Max - Min))
            {
                return Max - Math.Sqrt((1.0d - probability) * (Max - Min) * (Max - MostLikely));
            }
            else if (Max - Min == 0d)
            {
                return MostLikely;
            }
            return double.NaN;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new Triangular(Min, Mode, Max);
        }


    }
}