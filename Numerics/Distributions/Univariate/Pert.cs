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

using Numerics.Data.Statistics;
using Numerics.Mathematics.Optimization;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Distributions
{

    /// <summary>
    /// The PERT probability distribution.
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
    ///     In probability and statistics, the PERT distribution is a family of continuous probability distributions
    ///     defined by the minimum (a), most likely (b) and maximum (c) values that a variable can take.
    ///     It is a transformation of the four-parameter Beta distribution.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// Wikipedia contributors, "PERT distribution,". Wikipedia, The Free
    /// Encyclopedia. Available at: <see href="https://en.wikipedia.org/wiki/PERT_distribution"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class Pert : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IBootstrappable
    {

        /// <summary>
        /// Constructs a PERT distribution with min = 0.0, max = 1.0, and mode = 0.5.
        /// </summary>
        public Pert()
        {
            SetParameters(0.0d, 0.5d, 1.0d);
        }

        /// <summary>
        /// Constructs a PERT distribution with specified min, max, and mode.
        /// </summary>
        /// <param name="min">The minimum possible value of the distribution.</param>
        /// <param name="mode">The mode, or most likely, value of the distribution.</param>
        /// <param name="max">The maximum possible value of the distribution.</param>
        public Pert(double min, double mode, double max)
        {
            SetParameters(min, mode, max);
        }

        private GeneralizedBeta _beta = new GeneralizedBeta();
        private double _min;
        private double _max;
        private double _mode;

        /// <summary>
        /// Gets the underlying Generalized Beta distribution.
        /// </summary>
        public GeneralizedBeta Beta => _beta;

        /// <summary>
        /// Get and set the min of the distribution.
        /// </summary>
        public double Min
        {
            get { return _min; }
            set { SetParameters(value, MostLikely, Max); }
        }

        /// <summary>
        /// Get and set the mode of the distribution.
        /// </summary>
        public double MostLikely
        {
            get { return _mode; }
            set { SetParameters(Min, value, Max); }
        }

        /// <summary>
        /// Get and set the max of the distribution.
        /// </summary>
        public double Max
        {
            get { return _max; }
            set { SetParameters(Min, MostLikely, value); }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 3; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.Pert; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "PERT"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "PERT"; }
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
            get 
            {
                if (_min.AlmostEquals(_max) && _min.AlmostEquals(_mode)) return _min;
                return (Min + 4 * MostLikely + Max) / 6; 
            }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get 
            {
                if (_min.AlmostEquals(_max) && _min.AlmostEquals(_mode)) return _min;
                return (Min + 6 * MostLikely + Max) / 8; 
            }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get { return MostLikely; }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get { return Math.Sqrt((Mean - Min) * (Max - Mean) / 7); }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get { return _beta.Skewness; }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get { return _beta.Kurtosis; }
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
            if (estimationMethod == ParameterEstimationMethod.MethodOfPercentiles)
            {
                var data = sample.ToList();
                data.Sort();
                var five = Statistics.Percentile(data, 0.05, true);
                var fifty = Statistics.Percentile(data, 0.5, true);
                var ninetyFive = Statistics.Percentile(data, 0.95, true);
                var pertP = new PertPercentile(five, fifty, ninetyFive);
                pertP.SolveParameters();
                var pert = pertP.ToPert();
                SetParameters(pert.Min, pert.Mode, pert.Max);
            }
            else if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
            {
                var min = Tools.Min(sample);
                var max = Tools.Max(sample);
                var mean = Tools.Mean(sample);
                var mode = (mean * 6d - max - min) / 4d;
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
            var newDistribution = (Pert)Clone();
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
            // Return PERT
            if (_parametersValid == true) _beta = GeneralizedBeta.PERT(_min, _mode, _max);
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
                if (throwException) throw new ArgumentOutOfRangeException(nameof(Min), "The min cannot be greater than the max.");
                return new ArgumentOutOfRangeException(nameof(Min), "The min cannot be greater than the max.");
            }
            else if (double.IsNaN(mode) || double.IsInfinity(mode) || mode < min || mode > max)
            {
                if (throwException) throw new ArgumentOutOfRangeException(nameof(MostLikely), "The mode (most likely) must be between the min and max.");
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
            var mode = (mean * 6d - max - min) / 4d;
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
                var N = new Pert();
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
            if (_parametersValid == false) ValidateParameters(Min, Mode, Max, true);
            // These checks are done specifically for an application where a 
            // user inputs min = max = mode
            if (_min.AlmostEquals(_max) && _min.AlmostEquals(_mode)) return 0.0d;
            if (double.IsNaN(_mode)) return 0.0d;
            return _beta.PDF(x);
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameters(Min, Mode, Max, true);
            // These checks are done specifically for an application where a 
            // user inputs min = max = mode
            if (_min.AlmostEquals(_max) && _min.AlmostEquals(_mode)) return 1d;
            if (double.IsNaN(_mode)) return 1;
            return _beta.CDF(x);
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameters(Min, Mode, Max, true);
            // These checks are done specifically for an application where a 
            // user inputs min = max = mode
            if (_min.AlmostEquals(_max) && _min.AlmostEquals(_mode)) return Min;
            if (double.IsNaN(_mode)) return Min;
            return _beta.InverseCDF(probability);
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new Pert(Min, Mode, Max);
        }


    }
}