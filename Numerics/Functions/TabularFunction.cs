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

        /// <summary>
        /// The extrapolation policy applied to out-of-range lookups. Default = None, which
        /// reproduces the historical endpoint hold exactly.
        /// </summary>
        /// <remarks>
        /// One policy governs both lookup directions of the same extended boundary segments:
        /// <see cref="Function(double)"/> extends on the x-axis sides, and
        /// <see cref="InverseFunction(double)"/> extends on the y-lookup sides. Extension is
        /// linear in the configured transform space (see
        /// <see cref="OrderedPairedData.GetYFromX(double, Transform, Transform, ExtrapolationSides)"/>),
        /// and the <see cref="AllowNegativeYValues"/> floor still binds on extended forward lookups.
        /// </remarks>
        public ExtrapolationSides Extrapolation { get; set; } = ExtrapolationSides.None;

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
            double y = opd.GetYFromX(x, XTransform, YTransform, Extrapolation);
            y = AllowNegativeYValues == false && (double.IsNaN(y) || y < 0) ? 0 : y;
            return y;
        }

        /// <inheritdoc/>
        public double InverseFunction(double y)
        {
            // Validate parameters
            if (ParametersValid == false) ValidateParameters(new double[] { 0 }, true);
            y = AllowNegativeYValues == false && (double.IsNaN(y) || y < 0) ? 0 : y;
            return opd.GetXFromY(y, XTransform, YTransform, Extrapolation);
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
            // Conditional presence: the attribute is written only when non-default so that every
            // pre-existing serialized form remains byte-identical.
            if (Extrapolation != ExtrapolationSides.None)
            {
                result.SetAttributeValue(nameof(Extrapolation), Extrapolation.ToString());
            }
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
        /// <exception cref="ArgumentException">Thrown when the element carries no embedded uncertain ordered paired data, or a present attribute is malformed, non-finite, or undefined.</exception>
        public static TabularFunction FromXElement(XElement xElement)
        {
            if (xElement == null) throw new ArgumentNullException(nameof(xElement));
            var tableElement = xElement.Element("UncertainOrderedPairedData");
            if (tableElement == null)
                throw new ArgumentException("The serialized tabular function is missing its embedded UncertainOrderedPairedData.", nameof(xElement));

            var function = new TabularFunction(new UncertainOrderedPairedData(tableElement));
            var xTransformAttribute = xElement.Attribute(nameof(XTransform));
            if (xTransformAttribute != null)
            {
                if (!Enum.TryParse(xTransformAttribute.Value, out Transform xTransform)
                    || !Enum.IsDefined(typeof(Transform), xTransform))
                    throw new ArgumentException("The serialized X transform is invalid.", nameof(xElement));
                function.XTransform = xTransform;
            }
            var yTransformAttribute = xElement.Attribute(nameof(YTransform));
            if (yTransformAttribute != null)
            {
                if (!Enum.TryParse(yTransformAttribute.Value, out Transform yTransform)
                    || !Enum.IsDefined(typeof(Transform), yTransform))
                    throw new ArgumentException("The serialized Y transform is invalid.", nameof(xElement));
                function.YTransform = yTransform;
            }
            var extrapolationAttribute = xElement.Attribute(nameof(Extrapolation));
            if (extrapolationAttribute != null)
            {
                if (!Enum.TryParse(extrapolationAttribute.Value, out ExtrapolationSides extrapolation)
                    || !Enum.IsDefined(typeof(ExtrapolationSides), extrapolation))
                    throw new ArgumentException("The serialized extrapolation policy is invalid.", nameof(xElement));
                function.Extrapolation = extrapolation;
            }
            var allowNegativeAttribute = xElement.Attribute(nameof(AllowNegativeYValues));
            if (allowNegativeAttribute != null)
            {
                if (!bool.TryParse(allowNegativeAttribute.Value, out bool allowNegative))
                    throw new ArgumentException("The serialized negative-Y policy is invalid.", nameof(xElement));
                function.AllowNegativeYValues = allowNegative;
            }
            if (TryReadFiniteDouble(xElement, nameof(Minimum), out double minimum))
                function.Minimum = minimum;
            if (TryReadFiniteDouble(xElement, nameof(Maximum), out double maximum))
                function.Maximum = maximum;
            return function;
        }

        /// <summary>Reads an optional finite double attribute.</summary>
        private static bool TryReadFiniteDouble(XElement xElement, string attributeName, out double value)
        {
            var attribute = xElement.Attribute(attributeName);
            if (attribute == null)
            {
                value = 0d;
                return false;
            }
            if (!double.TryParse(attribute.Value, NumberStyles.Any, CultureInfo.InvariantCulture, out value)
                || !Tools.IsFinite(value))
            {
                throw new ArgumentException("The serialized " + attributeName + " value is invalid.", nameof(xElement));
            }
            return true;
        }
    }
}
