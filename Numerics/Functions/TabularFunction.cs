using Numerics.Data;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Xml.Linq;

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
        public ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            var errors = PairedData.GetErrors();
            if (errors.Count > 0)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(PairedData), "The uncertain ordered paired data has errors.");
                return new ArgumentOutOfRangeException(nameof(PairedData), "The uncertain ordered paired data has errors.");
            }
            return null!;
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

        /// <summary>
        /// Serializes the function's configuration to an XElement: the embedded uncertain ordered
        /// paired data plus the axis transforms, the negative-Y policy, and the support bounds.
        /// <see cref="ConfidenceLevel"/> is runtime sampling state and is deliberately not
        /// serialized; <see cref="IsDeterministic"/> derives from the table's distribution type.
        /// </summary>
        /// <returns>An XElement representation of the tabular function.</returns>
        public XElement ToXElement()
        {
            var result = new XElement(nameof(TabularFunction));
            result.SetAttributeValue(nameof(XTransform), XTransform.ToString());
            result.SetAttributeValue(nameof(YTransform), YTransform.ToString());
            result.SetAttributeValue(nameof(AllowNegativeYValues), AllowNegativeYValues.ToString());
            result.SetAttributeValue(nameof(Minimum), Minimum.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(Maximum), Maximum.ToString("G17", CultureInfo.InvariantCulture));
            result.Add(PairedData.SaveToXElement());
            return result;
        }

        /// <summary>
        /// Deserializes a tabular function from an XElement produced by <see cref="ToXElement"/>.
        /// </summary>
        /// <param name="xElement">The XElement to deserialize.</param>
        /// <returns>A new <see cref="TabularFunction"/>.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="xElement"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the element carries no embedded uncertain ordered paired data.</exception>
        public static TabularFunction FromXElement(XElement xElement)
        {
            if (xElement == null) throw new ArgumentNullException(nameof(xElement));
            var tableElement = xElement.Element("UncertainOrderedPairedData");
            if (tableElement == null)
                throw new ArgumentException("The serialized tabular function is missing its embedded UncertainOrderedPairedData.", nameof(xElement));

            var function = new TabularFunction(new UncertainOrderedPairedData(tableElement));
            if (Enum.TryParse(xElement.Attribute(nameof(XTransform))?.Value, out Transform xTransform))
                function.XTransform = xTransform;
            if (Enum.TryParse(xElement.Attribute(nameof(YTransform))?.Value, out Transform yTransform))
                function.YTransform = yTransform;
            if (bool.TryParse(xElement.Attribute(nameof(AllowNegativeYValues))?.Value, out bool allowNegative))
                function.AllowNegativeYValues = allowNegative;
            if (double.TryParse(xElement.Attribute(nameof(Minimum))?.Value, NumberStyles.Any, CultureInfo.InvariantCulture, out double minimum))
                function.Minimum = minimum;
            if (double.TryParse(xElement.Attribute(nameof(Maximum))?.Value, NumberStyles.Any, CultureInfo.InvariantCulture, out double maximum))
                function.Maximum = maximum;
            return function;
        }
    }
}
