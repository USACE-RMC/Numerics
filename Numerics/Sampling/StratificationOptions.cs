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
using System.Xml.Linq;

namespace Numerics.Sampling
{

    /// <summary>
    /// A class for stratification options.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class StratificationOptions
    {

        /// <summary>
        /// Get the lower bound of the starting bin.
        /// </summary>
        public double LowerBound { get; private set; }

        /// <summary>
        /// Get the upper bound of the last bin.
        /// </summary>
        public double UpperBound { get; private set; }

        /// <summary>
        /// Returns the number of bins. Must be greater than 1.
        /// </summary>
        public int NumberOfBins { get; private set; }

        /// <summary>
        /// Determines if the values to be stratified are probabilities [0,1].
        /// </summary>
        public bool IsProbability { get; private set; }

        /// <summary>
        /// Determines if the stratification options are valid.
        /// </summary>
        public bool IsValid { get; private set; }

        /// <summary>
        /// Constructs a stratification options class.
        /// </summary>
        /// <param name="lowerBound">The lower bound of the starting stratification bin.</param>
        /// <param name="upperBound">The upper bound of the last stratification bin.</param>
        /// <param name="numberOfBins">The number of stratification bins.</param>
        /// <param name="isProbability">Optional. Determines if the values to be stratified are probabilities [0,1]. Default = False.</param>
        public StratificationOptions(double lowerBound, double upperBound, int numberOfBins, bool isProbability = false)
        {
            // Set the options
            LowerBound = lowerBound;
            UpperBound = upperBound;
            NumberOfBins = numberOfBins;
            IsProbability = isProbability;
            Validate();
        }

        /// <summary>
        /// Initialize a new instance of the stratification options class from XElement.
        /// </summary>
        /// <param name="element">XElement to deserialize into a stratification options class.</param>
        public StratificationOptions(XElement element)
        {
            // Get required data
            var lowerBoundAttr = element.Attribute(nameof(LowerBound));
            if (lowerBoundAttr != null)
            {
                double.TryParse(lowerBoundAttr.Value, NumberStyles.Any, CultureInfo.InvariantCulture, out var lower);
                LowerBound = lower;
            }

            var upperBoundAttr = element.Attribute(nameof(UpperBound));
            if (upperBoundAttr != null)
            {
                double.TryParse(upperBoundAttr.Value, NumberStyles.Any, CultureInfo.InvariantCulture, out var upper);
                UpperBound = upper;
            }

            var numberOfBinsAttr = element.Attribute(nameof(NumberOfBins));
            if (numberOfBinsAttr != null)
            {
                int.TryParse(numberOfBinsAttr.Value, NumberStyles.Any, CultureInfo.InvariantCulture, out var nBins);
                NumberOfBins = nBins;
            }

            var isProbabilityAttr = element.Attribute(nameof(IsProbability));
            if (isProbabilityAttr != null)
            {
                bool.TryParse(isProbabilityAttr.Value, out var isProbability);
                IsProbability = isProbability;
            }

            Validate();
        }

        /// <summary>
        /// Returns an XElement of stratification options, can be used for serialization.
        /// </summary>
        public XElement SaveToXElement()
        {
            var result = new XElement("StratificationOptions");
            result.SetAttributeValue(nameof(LowerBound), LowerBound.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(UpperBound), UpperBound.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(NumberOfBins), NumberOfBins.ToString(CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(IsProbability), IsProbability.ToString());
            return result;
        }

        /// <summary>
        /// Returns an XElement of a list of stratification options, can be used for serialization.
        /// </summary>
        /// <param name="stratificationOptions">Collection of stratification options.</param>
        public static XElement SaveToXElement(IList<StratificationOptions> stratificationOptions)
        {
            var result = new XElement("StratificationOptionsCollection");
            if (stratificationOptions == null) return result;
            foreach (var s in stratificationOptions)
                result.Add(s.SaveToXElement());
            return result;
        }

        /// <summary>
        /// Converts an XElement to a list of stratification options.
        /// </summary>
        /// <param name="element">the XElement that will be deserialized to a list of stratification options.</param>
        public static List<StratificationOptions> XElementToStratificationOptions(XElement element)
        {
            var result = new List<StratificationOptions>();
            if (element == null) return result;
            foreach (XElement s in element.Elements("StratificationOptions"))
                result.Add(new StratificationOptions(s));
            return result;
        }

        /// <summary>
        /// Deep clones the stratification options object.
        /// </summary>
        /// <returns>A deep clone of the stratification options object.</returns>
        public StratificationOptions Clone()
        {
            return new StratificationOptions(LowerBound, UpperBound, NumberOfBins, IsProbability);
        }

        /// <summary>
        /// Deep clones a list of stratification options to a list of stratification options.
        /// </summary>
        /// <param name="stratificationOptions">Collection of stratification options to be deep cloned.</param>
        public static List<StratificationOptions> Clone(List<StratificationOptions> stratificationOptions)
        {
            var result = new List<StratificationOptions>();
            if (stratificationOptions == null) return result;
            foreach (var valueInCollection in stratificationOptions)
                result.Add(valueInCollection.Clone());
            return result;
        }

        /// <summary>
        /// Gets a list of errors for the stratification options.
        /// </summary>
        public List<string> GetErrors()
        {
            var result = new List<string>();
            if (IsValid == true)
                return result;
            if (LowerBound >= UpperBound)
            {
                result.Add("The upper bound '" + UpperBound + "' must be greater than the lower bound '" + LowerBound + "'.");
            }

            if (NumberOfBins < 2.0d)
            {
                result.Add("The number of bins '" + NumberOfBins + "' must be greater than 1.");
            }

            if (IsProbability == true && LowerBound < 0.0d)
            {
                result.Add("The lower bound must be greater than or equal to 0.");
            }

            if (IsProbability == true && UpperBound > 1.0d)
            {
                result.Add("The upper bound must be less than or equal to 1.");
            }

            return result;
        }

        /// <summary>
        /// Validates the options and sets the IsValid flag.
        /// </summary>
        private void Validate()
        {
            IsValid = GetErrors().Count == 0;
        }

        /// <summary>
        /// Compares objects for equality.
        /// </summary>
        public override bool Equals(object? obj)
        {
            if (obj is not StratificationOptions other)
                return false;

            return LowerBound == other.LowerBound &&
                   UpperBound == other.UpperBound &&
                   NumberOfBins == other.NumberOfBins &&
                   IsProbability == other.IsProbability;
        }

        /// <summary>
        /// Serves as the default hash function.
        /// </summary>
        /// <returns>A hash code for the current object.</returns>
        public override int GetHashCode()
        {
            return HashCode.Combine(LowerBound, UpperBound, NumberOfBins, IsProbability);
        }
    }
}