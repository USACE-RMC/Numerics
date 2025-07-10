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
using System.Globalization;
using System.Reflection;
using System.Xml.Linq;
using Numerics.Mathematics.Optimization;

namespace Numerics.Distributions
{
    /// <summary>
    /// The Pert percentile distribution.
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
    ///     This version of the PERT is parameterized using the 5th, 50th, and 95th percentiles.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/PERT_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class PertPercentile : UnivariateDistributionBase
    {

        /// <summary>
        /// Constructs a PERT distribution with 5th = 0.05, 50th = 0.5, and 95th = 0.95.
        /// </summary>
        public PertPercentile()
        {
            SetParameters(0.05d, 0.5d, 0.95d);
        }

        /// <summary>
        /// Constructs a PERT distribution with specified 5th, 50th, and 95th percentiles.
        /// </summary>
        /// <param name="fifth">The 5th percentile value of the distribution.</param>
        /// <param name="fiftieth">The 50th percentile value of the distribution.</param>
        /// <param name="ninetyFifth">The 95th percentile value of the distribution.</param>
        public PertPercentile(double fifth, double fiftieth, double ninetyFifth)
        {
            SetParameters(fifth, fiftieth, ninetyFifth);
        }

        private bool _parametersSolved = false;
        private GeneralizedBeta _beta = new GeneralizedBeta();
        private double _5th;
        private double _50th;
        private double _95th;

        /// <summary>
        /// Gets and sets the 5th percentile.
        /// </summary>
        public double Percentile5th
        {
            get { return _5th; }
            set
            {
                if (_5th != value)
                {
                    _5th = value;
                    // validate parameters
                    _parametersValid = ValidateParameters(_5th, _50th, _95th, false) is null;
                    _parametersSolved = false;
                }
            }
        }

        /// <summary>
        /// Gets and sets the 50th percentile.
        /// </summary>
        public double Percentile50th
        {
            get { return _50th; }
            set
            {
                if (_50th != value)
                {
                    _50th = value;
                    // validate parameters
                    _parametersValid = ValidateParameters(_5th, _50th, _95th, false) is null;
                    _parametersSolved = false;
                }
            }
        }

        /// <summary>
        /// Gets and sets the 95th percentile.
        /// </summary>
        public double Percentile95th
        {
            get { return _95th; }
            set
            {
                if (_95th != value)
                {
                    _95th = value;
                    // validate parameters
                    _parametersValid = ValidateParameters(_5th, _50th, _95th, false) is null;
                    _parametersSolved = false;
                }
            }
        }

        /// <summary>
        /// The minimum allowable value that can be sampled.
        /// </summary>
        public double MinAllowableValue { get; set; } = double.NegativeInfinity;

        /// <summary>
        /// The maximum allowable value that can be sampled.
        /// </summary>
        public double MaxAllowableValue { get; set; } = double.PositiveInfinity;


        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 3; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.PertPercentile; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "PERT-Percentile"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "PERT-%"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[3, 2];
                parmString[0, 0] = "5%";
                parmString[1, 0] = "50%";
                parmString[2, 0] = "95%";
                parmString[0, 1] = Percentile5th.ToString();
                parmString[1, 1] = Percentile50th.ToString();
                parmString[2, 1] = Percentile95th.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["5%", "50%", "95%"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Percentile5th), nameof(Percentile50th), nameof(Percentile95th)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Percentile5th, Percentile50th, Percentile95th]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get 
            {
                if (_beta.Min.AlmostEquals(_beta.Max)) return _beta.Min;
                return _beta.Mean < MinAllowableValue ? MinAllowableValue : _beta.Mean > MaxAllowableValue ? MaxAllowableValue : _beta.Mean;  
            }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get 
            {
                if (_beta.Min.AlmostEquals(_beta.Max)) return _beta.Min;
                return _beta.Median < MinAllowableValue ? MinAllowableValue : _beta.Median > MaxAllowableValue ? MaxAllowableValue : _beta.Median; 
            }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get 
            {
                if (_beta.Min.AlmostEquals(_beta.Max)) return _beta.Min;
                return _beta.Mode < MinAllowableValue ? MinAllowableValue : _beta.Mode > MaxAllowableValue ? MaxAllowableValue : _beta.Mode; 
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get { return _beta.StandardDeviation; }
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
            get { return Math.Max(_beta.Minimum, MinAllowableValue); }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return Math.Min(_beta.Maximum, MaxAllowableValue); }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [double.NegativeInfinity, Percentile5th, Percentile50th]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [Percentile50th, Percentile95th, double.PositiveInfinity]; }
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="fifth">The 5th percentile value of the distribution.</param>
        /// <param name="fiftieth">The 50th percentile value of the distribution.</param>
        /// <param name="ninetyFifth">The 95th percentile value of the distribution.</param>
        public void SetParameters(double fifth, double fiftieth, double ninetyFifth)
        {
            Percentile5th = fifth;
            Percentile50th = fiftieth;
            Percentile95th = ninetyFifth;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1], parameters[2]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="fifth">The 5th percentile value of the distribution.</param>
        /// <param name="fiftieth">The 50th percentile value of the distribution.</param>
        /// <param name="ninetyFifth">The 95th percentile value of the distribution.</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        private ArgumentOutOfRangeException ValidateParameters(double fifth, double fiftieth, double ninetyFifth, bool throwException)
        {
            if (double.IsNaN(fifth) || double.IsInfinity(fifth) ||
                double.IsNaN(ninetyFifth) || double.IsInfinity(ninetyFifth) || fifth > ninetyFifth)
            {
                if (throwException) throw new ArgumentOutOfRangeException(nameof(Percentile5th), "The 5% cannot be greater than the 95%.");
                return new ArgumentOutOfRangeException(nameof(Percentile5th), "The 5% cannot be greater than the 95%.");
            }
            else if (double.IsNaN(fiftieth) || double.IsInfinity(fiftieth) || fiftieth < fifth || fiftieth > ninetyFifth)
            {
                if (throwException) throw new ArgumentOutOfRangeException(nameof(Percentile50th), "The 50% must be between the 5% and 95%.");
                return new ArgumentOutOfRangeException(nameof(Percentile50th), "The 50% must be between the 5% and 95%.");
            }
            return null;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            return ValidateParameters(parameters[0], parameters[1], parameters[2], throwException);
        }


        /// <summary>
        /// Solve the PERT parameters (min, mode, max) given the percentiles.
        /// </summary>
        public void SolveParameters()
        {
            if (_parametersValid == false) return;
            if (_parametersSolved == true) return;

            double fifth = Percentile5th;
            double fiftieth = Percentile50th;
            double ninetyFifth = Percentile95th;

            if (fifth == fiftieth && fiftieth == ninetyFifth)
            {
                _beta = GeneralizedBeta.PERT(fifth, fiftieth, ninetyFifth);
                _parametersSolved = true;
                return;
            }

            _beta = GeneralizedBeta.PERT(fifth, fiftieth, ninetyFifth);
            // Set wide bounds
            double min = fifth - (ninetyFifth - fifth) * 2;
            double max = ninetyFifth + (ninetyFifth - fifth) * 2;
            var Initials = new double[] { _beta.Min, _beta.Mode, _beta.Max };
            var Lowers = new double[] { min, min, min };
            var Uppers = new double[] { max, max, max };

            // Objective function is the sum of squared errors
            double sse(double[] x)
            {
                var dist = new Pert();
                try
                {
                    dist = new Pert(x[0], x[1], x[2]);
                    if (dist.ParametersValid == false) return double.MaxValue;

                    double SSE = 0d;
                    SSE += Math.Pow(fifth - dist.InverseCDF(0.05), 2d);
                    SSE += Math.Pow(fiftieth - dist.InverseCDF(0.5), 2d);
                    SSE += Math.Pow(ninetyFifth - dist.InverseCDF(0.95), 2d);
                    return SSE;
                }
                catch
                {
                    return double.MaxValue;
                }

            }

            var solver = new NelderMead(sse, 3, Initials, Lowers, Uppers);
            solver.MaxFunctionEvaluations = 200;
            solver.RelativeTolerance = 1E-8;
            solver.AbsoluteTolerance = 1E-8;
            solver.ComputeHessian = false;
            solver.ReportFailure = false;
            solver.Minimize();
            var solution = solver.BestParameterSet.Values;
            _beta = GeneralizedBeta.PERT(solution[0], solution[1], solution[2]);
            _parametersSolved = true;
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            if (_parametersValid == false) ValidateParameters(Percentile5th, Percentile50th, Percentile95th, true);
            if (_parametersSolved == false) SolveParameters();
            // These checks are done specifically for an application where a 
            // user inputs min = max = mode
            if (_beta.Min.AlmostEquals(_beta.Max) && _beta.Min.AlmostEquals(_beta.Mode)) return 0.0d;
            if (double.IsNaN(_beta.Mode)) return 0.0d;
            //
            if (x < MinAllowableValue) x = MinAllowableValue;
            if (x > MaxAllowableValue) x = MaxAllowableValue;
            return _beta.PDF(x);
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            if (_parametersValid == false) ValidateParameters(Percentile5th, Percentile50th, Percentile95th, true);
            if (_parametersSolved == false) SolveParameters();
            // These checks are done specifically for an application where a 
            // user inputs min = max = mode
            if (_beta.Min.AlmostEquals(_beta.Max) && _beta.Min.AlmostEquals(_beta.Mode)) return 1d;
            if (double.IsNaN(_beta.Mode)) return 1d;
            //
            if (x < MinAllowableValue) x = MinAllowableValue;
            if (x > MaxAllowableValue) x = MaxAllowableValue;
            return _beta.CDF(x);
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameters(Percentile5th, Percentile50th, Percentile95th, true);
            if (_parametersSolved == false) SolveParameters();
            //
            var x = _beta.Min;
            if (x < MinAllowableValue) x = MinAllowableValue;
            if (x > MaxAllowableValue) x = MaxAllowableValue;
            // These checks are done specifically for an application where a 
            // user inputs min = max = mode
            if (_beta.Min.AlmostEquals(_beta.Max) && _beta.Min.AlmostEquals(_beta.Mode)) return x;
            if (double.IsNaN(_beta.Mode)) return x;
            //
            try
            {
                x = _beta.InverseCDF(probability);
            }
            catch (ArithmeticException)
            {
                return x;
            } 
            //
            if (x < MinAllowableValue) x = MinAllowableValue;
            if (x > MaxAllowableValue) x = MaxAllowableValue;
            return x;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new PertPercentile(Percentile5th, Percentile50th, Percentile95th) { MinAllowableValue = MinAllowableValue,
                MaxAllowableValue = MaxAllowableValue,
                _parametersSolved = _parametersSolved,
                _parametersValid = _parametersValid,
                _beta = (GeneralizedBeta)_beta.Clone()};
        }

        /// <summary>
        /// Return the Pert-Percentile as a Pert distribution.
        /// </summary>
        public Pert ToPert()
        {
            return new Pert(_beta.Min, _beta.Mode, _beta.Max);
        }

        /// <inheritdoc/>
        public override XElement ToXElement()
        {
            var result = new XElement("Distribution");
            result.SetAttributeValue(nameof(Type), Type.ToString());
            result.SetAttributeValue("ParametersSolved", _parametersSolved.ToString());
            // Parameters
            var parms = GetParameters;
            var parmStrings = new string[NumberOfParameters];
            for (int i = 0; i < NumberOfParameters; i++)
            {
                parmStrings[i] = parms[i].ToString("G17", CultureInfo.InvariantCulture);
            }
            result.SetAttributeValue("Parameters", String.Join("|", parmStrings));
            // Beta Distribution parameters
            parms = _beta.GetParameters;
            parmStrings = new string[_beta.NumberOfParameters];
            for (int i = 0; i < _beta.NumberOfParameters; i++)
            {
                parmStrings[i] = parms[i].ToString("G17", CultureInfo.InvariantCulture);
            }
            result.SetAttributeValue("BetaParameters", String.Join("|", parmStrings));
            return result;
        }

        /// <summary>
        /// Create a Pert-Percentile distribution from XElement.
        /// </summary>
        /// <param name="xElement">The XElement to deserialize.</param>
        /// <returns>A new mixture distribution.</returns>
        public static PertPercentile FromXElement(XElement xElement)
        {
            UnivariateDistributionType type = UnivariateDistributionType.Deterministic;
            if (xElement.Attribute(nameof(UnivariateDistributionBase.Type)) != null)
            {
                Enum.TryParse(xElement.Attribute(nameof(UnivariateDistributionBase.Type)).Value, out type);

            }
            if (type == UnivariateDistributionType.PertPercentile)
            {
                bool parametersSolved = false;
                if (xElement.Attribute("ParametersSolved") != null)
                {
                    bool.TryParse(xElement.Attribute("ParametersSolved").Value, out parametersSolved);
                }
                else
                {
                    // This was saved from Base method
                    var dist = new PertPercentile();
                    var names = dist.GetParameterPropertyNames;
                    var parms = dist.GetParameters;
                    var vals = new double[dist.NumberOfParameters];
                    for (int i = 0; i < dist.NumberOfParameters; i++)
                    {
                        if (xElement.Attribute(names[i]) != null)
                        {
                            double.TryParse(xElement.Attribute(names[i]).Value, NumberStyles.Any, CultureInfo.InvariantCulture, out vals[i]);
                        }
                    }
                    dist.SetParameters(vals);
                    return dist;

                }

                var beta = new GeneralizedBeta();
                if (xElement.Attribute("BetaParameters") != null)
                {
                    var vals = xElement.Attribute("BetaParameters").Value.Split('|');
                    var parameters = new List<double>();
                    for (int i = 0; i < vals.Length; i++)
                    {
                        double.TryParse(vals[i], NumberStyles.Any, CultureInfo.InvariantCulture, out var parm);
                        parameters.Add(parm);
                    }
                    beta.SetParameters(parameters);
                }
                double _5th = 0, _50th = 0, _95th = 0;
                if (xElement.Attribute("Parameters") != null)
                {
                    var vals = xElement.Attribute("Parameters").Value.Split('|');
                    double.TryParse(vals[0], NumberStyles.Any, CultureInfo.InvariantCulture, out _5th);
                    double.TryParse(vals[1], NumberStyles.Any, CultureInfo.InvariantCulture, out _50th);
                    double.TryParse(vals[2], NumberStyles.Any, CultureInfo.InvariantCulture, out _95th);
                }

                var pert = new PertPercentile() 
                { 
                    _parametersSolved = parametersSolved, 
                    _beta = beta, 
                    _5th = _5th, 
                    _50th = _50th, 
                    _95th = _95th };

                return pert;
            }
            else
            {
                return null;
            }
        }

    }
}
