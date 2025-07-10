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

using Numerics.Data;
using System;
using System.Collections.Generic;

namespace Numerics.Functions
{
    /// <summary>
    /// A class for a tabular, or nonparametric, function.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class TabularFunction : IUnivariateFunction
    {

        /// <summary>
        /// Constructs a tabular function from uncertain ordered paired data.
        /// </summary>
        /// <param name="pairedData">The uncertain ordered paired data.</param>
        public TabularFunction(UncertainOrderedPairedData pairedData)
        {
            _pairedData = pairedData;
            opd = _pairedData.CurveSample();
        }

        private UncertainOrderedPairedData _pairedData;
        private OrderedPairedData opd;
        private double _confidenceLevel = -1;

        /// <summary>
        /// The uncertain ordered paired data. 
        /// </summary>
        public UncertainOrderedPairedData PairedData => _pairedData;

        /// <summary>
        /// The transform for the x-values. Default = None.
        /// </summary>
        public Transform XTransform { get; set; } = Transform.None;

        /// <summary>
        /// The transform for the y-values. Default = None.
        /// </summary>
        public Transform YTransform { get; set; } = Transform.None;

        /// <inheritdoc/>
        public int NumberOfParameters => 1;

        /// <inheritdoc/>
        public bool ParametersValid => PairedData.IsValid;

        /// <inheritdoc/>
        public double Minimum { get; set; } = double.MinValue;

        /// <inheritdoc/>
        public double Maximum { get; set; } = double.MaxValue;

        /// <inheritdoc/>
        public double[] MinimumOfParameters => new double[] { double.MinValue };

        /// <inheritdoc/>
        public double[] MaximumOfParameters => new double[] { double.MaxValue};

        /// <inheritdoc/>
        public bool IsDeterministic 
        {
            get
            { 
                if (_pairedData.Distribution == Distributions.UnivariateDistributionType.Deterministic)
                {
                    return true;
                }
                else
                {
                    return false;
                }                
            }
            set
            {
                if (value == true)
                {
                    _pairedData = new UncertainOrderedPairedData(_pairedData, _pairedData.StrictX, _pairedData.OrderX, _pairedData.StrictY, _pairedData.OrderY, Distributions.UnivariateDistributionType.Deterministic);
                }
            }
        }

        /// <inheritdoc/>
        public double ConfidenceLevel 
        { 
            get { return _confidenceLevel; }
            set
            {
                _confidenceLevel = value;
                if (_confidenceLevel < 0 )
                {
                    opd = _pairedData.CurveSample();
                }
                else
                {
                    opd = _pairedData.CurveSample(_confidenceLevel);
                }
                
            }
        }

        /// <summary>
        /// Determines if the tabular function can return negative Y values.
        /// </summary>
        public bool AllowNegativeYValues { get; set; } = true;

        /// <inheritdoc/>
        public void SetParameters(IList<double> parameters)
        {
            // This method is not implemented since the tabular function uses
            // uncertain paired data as the input. 
            throw new NotImplementedException();
        }

        /// <inheritdoc/>
        public ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            var errors = PairedData.GetErrors();
            if (errors.Count > 0)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(PairedData), "The uncertain ordered paired data has errors.");
                return new ArgumentOutOfRangeException(nameof(PairedData), "The uncertain ordered paired data has errors.");
            }
            return null;
        }

        /// <inheritdoc/>
        public double Function(double x)
        {
            // Validate parameters
            if (ParametersValid == false) ValidateParameters(new double[] {0}, true);
            double y = opd.GetYFromX(x, XTransform, YTransform);
            y = AllowNegativeYValues == false && (double.IsNaN(y) || y < 0) ? 0 : y;
            return y;
        }

        /// <inheritdoc/>
        public double InverseFunction(double y)
        {
            // Validate parameters
            if (ParametersValid == false) ValidateParameters(new double[] { 0 }, true);
            y = AllowNegativeYValues == false && (double.IsNaN(y) || y < 0) ? 0 : y;          
            return opd.GetXFromY(y, XTransform, YTransform);
        }
    }
}
