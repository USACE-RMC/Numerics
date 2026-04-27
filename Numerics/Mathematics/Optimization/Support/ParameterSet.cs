using System;
using System.Globalization;
using System.Linq;
using System.Xml.Linq;

namespace Numerics.Mathematics.Optimization
{
    /// <summary>
    /// A class for storing an optimization trial parameter set.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public struct ParameterSet
    {
        /// <summary>
        /// Constructs an empty parameter set.
        /// </summary>
        public ParameterSet() { }

        /// <summary>
        /// Constructs a parameter set.
        /// </summary>
        /// <param name="values">The parameter values.</param>
        /// <param name="fitness">The objective function result (or fitness) given the parameter set.</param>
        public ParameterSet(double[] values, double fitness)
        {
            Values = values;
            Fitness = fitness;
        }

        /// <summary>
        /// Constructs a parameter set.
        /// </summary>
        /// <param name="values">The parameter values.</param>
        /// <param name="fitness">The objective function result (or fitness) given the parameter set.</param>
        /// <param name="weight">The weight given to the parameter set.</param>
        public ParameterSet(double[] values, double fitness, double weight)
        {
            Values = values;
            Fitness = fitness;
            Weight = weight;
        }

        /// <summary>
        /// Constructs a parameter set.
        /// </summary>
        /// <param name="xElement">The XElement to deserialize.</param>
        public ParameterSet(XElement xElement)
        {

            var valuesAttr = xElement.Attribute(nameof(Values));
            if (valuesAttr != null && !string.IsNullOrWhiteSpace(valuesAttr.Value))
            {        
                var vals = valuesAttr.Value.Split('|');
                Values = new double[vals.Length];
                for (int i = 0; i < vals.Length; i++)
                {
                    double.TryParse(vals[i], NumberStyles.Any, CultureInfo.InvariantCulture, out var outVal);
                    Values[i] = outVal;
                }
            }
            var fitnessAttr = xElement.Attribute(nameof(Fitness));
            if (fitnessAttr != null)
            {
                double.TryParse(fitnessAttr.Value, NumberStyles.Any, CultureInfo.InvariantCulture, out var fitness);
                Fitness = fitness;
            }
            var weightAttr = xElement.Attribute(nameof(Weight));
            if (weightAttr != null)
            {
                double.TryParse(weightAttr.Value, NumberStyles.Any, CultureInfo.InvariantCulture, out var weight);
                Weight = weight;
            }
        }

        /// <summary>
        /// The trial parameter set values.
        /// </summary>
        public double[] Values = null!;

        /// <summary>
        /// The objective function result (or fitness) given the trial parameter set. 
        /// </summary>
        public double Fitness;

        /// <summary>
        /// An optional weight given to the parameter set values.
        /// </summary>
        public double Weight;

        /// <summary>
        /// Returns a clone of the point.
        /// </summary>
        /// <param name="deep">True to use a deep clone of the values array; false to use the array reference directly. Default = true.</param>
        public ParameterSet Clone(bool deep = true)
        {
            return deep ? new ParameterSet((double[])Values.Clone(), Fitness, Weight) : new ParameterSet(Values, Fitness, Weight);
        }

        /// <summary>
        /// Returns the matrix as XElement.
        /// </summary>
        public XElement ToXElement()
        {
            var result = new XElement(nameof(ParameterSet));
            if (Values != null)
                result.SetAttributeValue(nameof(Values), string.Join("|", Values.Select(v => v.ToString("G17", CultureInfo.InvariantCulture))));
            result.SetAttributeValue(nameof(Fitness), Fitness.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(Weight), Weight.ToString("G17", CultureInfo.InvariantCulture));
            return result;
        }
    }
}
