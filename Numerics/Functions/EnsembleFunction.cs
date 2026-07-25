using Numerics.Mathematics.Optimization;
using System;
using System.Collections.Generic;
using System.Xml.Linq;

namespace Numerics.Functions
{
    /// <summary>
    /// A posterior ensemble of a univariate function: a template function plus an array of
    /// posterior parameter sets, sampled by index or percentile as fresh configured clones.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// This is how an imported fitted function (e.g., a BestFit rating curve) carries knowledge
    /// uncertainty into a simulation engine. Both sampling surfaces are <b>pure</b>: every call
    /// returns a brand-new clone of the template configured with the selected parameter set, so
    /// concurrent realizations never share mutable state — the thread-safe alternative to
    /// mutating one instance's <c>ConfidenceLevel</c>/<c>SetParameters</c> inside a
    /// <c>Parallel.For</c> hot loop. The template is captured as an immutable serialized
    /// snapshot at construction (each clone re-materializes through
    /// <see cref="UnivariateFunctionFactory.CreateFromXElement(XElement)"/>), so later edits to
    /// the original instance never leak into samples.
    /// </para>
    /// <para>
    /// Percentile sampling maps u ∈ [0, 1] onto the index ladder as
    /// min(⌊u·N⌋, N − 1) — the standard posterior-lookup convention.
    /// </para>
    /// </remarks>
    [Serializable]
    public class EnsembleFunction
    {

        /// <summary>
        /// Construct a new ensemble function from a template and its posterior parameter sets.
        /// </summary>
        /// <param name="template">The template function (a library function type).</param>
        /// <param name="parameterSets">The posterior parameter sets; each must carry one value per template parameter.</param>
        /// <exception cref="ArgumentNullException">Thrown when either argument is null.</exception>
        /// <exception cref="ArgumentException">Thrown when no parameter sets are supplied, or a set's length does not match the template.</exception>
        /// <exception cref="NotSupportedException">Thrown when the template is not a serializable library function type.</exception>
        public EnsembleFunction(IUnivariateFunction template, IList<ParameterSet> parameterSets)
        {
            if (template == null) throw new ArgumentNullException(nameof(template));
            if (parameterSets == null) throw new ArgumentNullException(nameof(parameterSets));
            if (parameterSets.Count == 0) throw new ArgumentException("At least one parameter set is required.", nameof(parameterSets));

            int expectedLength = template.NumberOfParameters;
            _parameterSets = new ParameterSet[parameterSets.Count];
            for (int i = 0; i < parameterSets.Count; i++)
            {
                if (parameterSets[i].Values == null || parameterSets[i].Values.Length != expectedLength)
                    throw new ArgumentException("Parameter set " + i + " must carry exactly " + expectedLength + " values (one per template parameter).", nameof(parameterSets));
                _parameterSets[i] = parameterSets[i];
            }

            // The immutable template snapshot: parsing a private string per sample guarantees
            // clones never share XML tree state across threads.
            _templateXml = SerializeTemplate(template).ToString(SaveOptions.DisableFormatting);
        }

        /// <summary>
        /// Construct a new ensemble function from a serialized template element and its
        /// posterior parameter sets (the deserialization path).
        /// </summary>
        /// <param name="templateElement">The template's serialized form.</param>
        /// <param name="parameterSets">The posterior parameter sets.</param>
        /// <exception cref="ArgumentNullException">Thrown when either argument is null.</exception>
        /// <exception cref="ArgumentException">Thrown when no parameter sets are supplied.</exception>
        private EnsembleFunction(XElement templateElement, IList<ParameterSet> parameterSets)
        {
            if (templateElement == null) throw new ArgumentNullException(nameof(templateElement));
            if (parameterSets == null) throw new ArgumentNullException(nameof(parameterSets));
            if (parameterSets.Count == 0) throw new ArgumentException("At least one parameter set is required.", nameof(parameterSets));
            _parameterSets = new ParameterSet[parameterSets.Count];
            for (int i = 0; i < parameterSets.Count; i++) _parameterSets[i] = parameterSets[i];
            _templateXml = templateElement.ToString(SaveOptions.DisableFormatting);
        }

        private readonly string _templateXml;
        private readonly ParameterSet[] _parameterSets;

        /// <summary>
        /// The number of posterior parameter sets.
        /// </summary>
        public int Count => _parameterSets.Length;

        /// <summary>
        /// The posterior parameter sets.
        /// </summary>
        public IReadOnlyList<ParameterSet> ParameterSets => _parameterSets;

        /// <summary>
        /// Samples the posterior by index: a fresh template clone configured with the indexed
        /// parameter set.
        /// </summary>
        /// <param name="index">The posterior index in [0, Count).</param>
        /// <returns>A new configured function instance.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="index"/> is outside [0, Count).</exception>
        public IUnivariateFunction Sample(int index)
        {
            if (index < 0 || index >= _parameterSets.Length)
                throw new ArgumentOutOfRangeException(nameof(index), "The posterior index must be within [0, Count).");
            var clone = UnivariateFunctionFactory.CreateFromXElement(XElement.Parse(_templateXml));
            clone.SetParameters(_parameterSets[index].Values);
            return clone;
        }

        /// <summary>
        /// Samples the posterior by percentile: u maps onto the index ladder as
        /// min(⌊u·N⌋, N − 1).
        /// </summary>
        /// <param name="percentile">The percentile u in [0, 1].</param>
        /// <returns>A new configured function instance.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="percentile"/> is outside [0, 1].</exception>
        public IUnivariateFunction Sample(double percentile)
        {
            if (percentile < 0d || percentile > 1d)
                throw new ArgumentOutOfRangeException(nameof(percentile), "The percentile must be between 0 and 1.");
            int index = (int)Math.Floor(percentile * _parameterSets.Length);
            if (index > _parameterSets.Length - 1) index = _parameterSets.Length - 1;
            return Sample(index);
        }

        /// <summary>
        /// Serializes the ensemble to an XElement: the template's serialized form plus every
        /// posterior parameter set.
        /// </summary>
        /// <returns>An XElement representation of the ensemble.</returns>
        public XElement ToXElement()
        {
            var result = new XElement(nameof(EnsembleFunction));
            var template = new XElement("Template");
            template.Add(XElement.Parse(_templateXml));
            result.Add(template);
            var sets = new XElement("ParameterSets");
            for (int i = 0; i < _parameterSets.Length; i++)
                sets.Add(_parameterSets[i].ToXElement());
            result.Add(sets);
            return result;
        }

        /// <summary>
        /// Deserializes an ensemble from an XElement produced by <see cref="ToXElement"/>.
        /// </summary>
        /// <param name="xElement">The XElement to deserialize.</param>
        /// <returns>A new <see cref="EnsembleFunction"/>.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="xElement"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the template or parameter sets are missing.</exception>
        public static EnsembleFunction FromXElement(XElement xElement)
        {
            if (xElement == null) throw new ArgumentNullException(nameof(xElement));
            var templateElement = xElement.Element("Template");
            XElement? template = null;
            if (templateElement != null)
            {
                foreach (var child in templateElement.Elements()) { template = child; break; }
            }
            if (template == null)
                throw new ArgumentException("The serialized ensemble function is missing its template.", nameof(xElement));

            var sets = new List<ParameterSet>();
            var setsElement = xElement.Element("ParameterSets");
            if (setsElement != null)
            {
                foreach (var child in setsElement.Elements())
                    sets.Add(new ParameterSet(child));
            }
            if (sets.Count == 0)
                throw new ArgumentException("The serialized ensemble function carries no parameter sets.", nameof(xElement));

            return new EnsembleFunction(template, sets);
        }

        /// <summary>
        /// Serializes a template through its concrete ToXElement method.
        /// </summary>
        /// <param name="function">The template function.</param>
        /// <returns>The template's serialized form.</returns>
        /// <exception cref="NotSupportedException">Thrown when the template is not a library function type.</exception>
        private static XElement SerializeTemplate(IUnivariateFunction function)
        {
            if (function is LinearFunction linear) return linear.ToXElement();
            if (function is PowerFunction power) return power.ToXElement();
            if (function is TabularFunction tabular) return tabular.ToXElement();
            if (function is SegmentedPowerFunction segmented) return segmented.ToXElement();
            if (function is CompositeFunction composite) return composite.ToXElement();
            throw new NotSupportedException("The template function type '" + function.GetType().Name + "' does not support serialization.");
        }

    }
}
