using System;
using System.Collections.Generic;
using System.Globalization;
using System.Xml.Linq;

namespace Numerics.Sampling
{

    /// <summary>
    /// A class for stratification bins.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Stratified_sampling" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class StratificationBin : IComparable<StratificationBin>, ICloneable
    {
    
        /// <summary>
        /// Construct new stratification bin.
        /// </summary>
        /// <param name="lowerBound">The lower bound of the bin.</param>
        /// <param name="upperBound">The upper bound of the bin.</param>
        /// <param name="weight">Optional. The weight or probability width of the bin. The weight does not have to be equal to the width.
        /// Default = -1, which will make the weight = width.</param>
        public StratificationBin(double lowerBound, double upperBound, double weight = -1)
        {
            // validate inputs
            if (lowerBound > upperBound)
            {
                throw new ArgumentOutOfRangeException(nameof(lowerBound), "The upper bound must be greater than or equal to the lower bound.");
            }

            LowerBound = lowerBound;
            UpperBound = upperBound;
            Weight = weight < 0 ? upperBound - lowerBound : weight;
        }


        /// <summary>
        /// Initialize a new instance of the stratified X values class.
        /// </summary>
        /// <param name="element">XElement to deserialized into a stratified x value class.</param>
        public StratificationBin(XElement element)
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

            var weightAttr = element.Attribute(nameof(Weight));
            if (weightAttr != null)
            {
                double.TryParse(weightAttr.Value, NumberStyles.Any, CultureInfo.InvariantCulture, out var weight);
                Weight = weight;
            }
        }
    
        /// <summary>
        /// Get the lower bound of the bin.
        /// </summary>
        public double LowerBound { get; private set; }

        /// <summary>
        /// Get the upper bound of the bin.
        /// </summary>
        public double UpperBound { get; private set; }

        /// <summary>
        /// Gets the midpoint of the bin.
        /// </summary>
        public double Midpoint => (UpperBound + LowerBound) / 2.0;

        /// <summary>
        /// The weight given to the stratification bin. This is often the same value as the bin width.
        /// However, end bins can be assigned different weights to ensure unity.
        /// </summary>
        public double Weight { get; set; }

        /// <summary>
        /// Checks if a value falls within the bin range.
        /// </summary>
        /// <param name="x">The value to test.</param>
        /// <returns>True if x is within the bin (inclusive lower, exclusive upper).</returns>
        public bool Contains(double x)
        {
            return x >= LowerBound && x < UpperBound;
        }

        /// <summary>
        /// Comparison of two bins. The bins cannot be overlapping.
        /// </summary>
        /// <param name="other">The bin to compare to.</param>
        /// <returns>
        /// 0 if the upper bound and lower bound are bit-for-bit equal.
        /// +1 if this bin is lower than the compared bin.
        /// -1 otherwise.
        /// </returns>
        public int CompareTo(StratificationBin? other)
        {
            if (other is null) return 1;
            if (UpperBound > other.LowerBound && LowerBound < other.UpperBound)
                throw new ArgumentException("The bins cannot be overlapping.", nameof(other));

            if (UpperBound.Equals(other.UpperBound) && LowerBound.Equals(other.LowerBound))
                return 0;

            return other.UpperBound <= LowerBound ? 1 : -1;
        }

        /// <summary>
        /// Creates a copy of the stratification bin.
        /// </summary>
        public object Clone()
        {
            return new StratificationBin(LowerBound, UpperBound, Weight);
        }

        /// <summary>
        /// Checks whether two stratification bins are equal.
        /// </summary>
        public override bool Equals(object? obj)
        {
            if (!(obj is StratificationBin))
                return false;
            StratificationBin bin = (StratificationBin)obj;
            return LowerBound.Equals(bin.LowerBound) && UpperBound.Equals(bin.UpperBound);
        }

        /// <summary>
        /// Serves as the default hash function. Purposefuly does not include the weight, because the 
        /// equals and compare to methods also do not.
        /// </summary>
        /// <returns>A hash code for the current object.</returns>
        public override int GetHashCode()
        {
            return HashCode.Combine(LowerBound.GetHashCode(), UpperBound.GetHashCode());
        }

        /// <summary>
        /// Returns an XElement of a stratification bin, can be used for serialization.
        /// </summary>
        public XElement SaveToXElement()
        {
            var result = new XElement("StratificationBin");
            result.SetAttributeValue(nameof(LowerBound), LowerBound.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(UpperBound), UpperBound.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(Weight), Weight.ToString("G17", CultureInfo.InvariantCulture));
            return result;
        }

        /// <summary>
        /// Returns an XElement of a list of stratification bins, can be used for serialization.
        /// </summary>
        /// <param name="stratificationBinList">Collection of stratification bins.</param>
        public static XElement SaveToXElement(IList<StratificationBin> stratificationBinList)
        {
            var result = new XElement("StratificationBinList");
            if (stratificationBinList == null) return result;
            foreach (var s in stratificationBinList)
                result.Add(s.SaveToXElement());
            return result;
        }

        /// <summary>
        /// Converts an XElement to a list of stratification bins.
        /// </summary>
        /// <param name="element">the XElement that will be deserialized to a list of stratification bins.</param>
        public static List<StratificationBin> XElementToStratificationBinList(XElement element)
        {
            var result = new List<StratificationBin>();
            if (element == null) return result;
            foreach (XElement s in element.Elements("StratificationBin"))
                result.Add(new StratificationBin(s));
            return result;
        }
    }

    
}