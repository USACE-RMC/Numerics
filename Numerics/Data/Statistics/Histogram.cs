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
using System.Linq;

namespace Numerics.Data.Statistics
{
    /// <summary>
    /// Create a histogram from a sample of data, which is a visual representation of the distribution of the data.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item>
    /// <see href = "https://en.wikipedia.org/wiki/Histogram" />
    /// </item>
    /// <item>
    /// This class is modeled after the histogram class in the Math.NET Numerics library: <see href="http://numerics.mathdotnet.com"/> 
    /// ></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class Histogram
    {
        /// <summary>
        /// A histogram bin.
        /// </summary>
        [Serializable]
        public class Bin : IComparable<Bin>, ICloneable
        {
            /// <summary>
            /// Initialize a new instance of the histogram bin class.
            /// </summary>
            /// <param name="lowerBound">The lower bound of the bin.</param>
            /// <param name="upperBound">The upper bound of the bin.</param>
            /// <param name="frequency">The frequency within the bin.</param>
            public Bin(double lowerBound, double upperBound, int frequency = 0)
            {
                // validate inputs
                if (lowerBound > upperBound)
                {
                    throw new ArgumentOutOfRangeException(nameof(LowerBound), "The upper bound must be greater than the lower bound.");
                }

                if (frequency < 0.0d)
                {
                    throw new ArgumentOutOfRangeException(nameof(Frequency), "The count must be non-negative.");
                }

                LowerBound = lowerBound;
                UpperBound = upperBound;
                Frequency = frequency;
            }

            /// <summary>
            /// Get and set the lower bound of the bin.
            /// </summary>
            public double LowerBound { get; set; }

            /// <summary>
            /// Get and set the upper bound of the bin.
            /// </summary>
            public double UpperBound { get; set; }

            /// <summary>
            /// Gets the midpoint of the bin.
            /// </summary>
            public double Midpoint
            {
                get { return (UpperBound + LowerBound) / 2d; }
            }

            /// <summary>
            /// Get and set the frequency of the bin.
            /// </summary>
            public int Frequency { get; set; }

            /// <summary>
            /// Comparison of two bins. The bins cannot be overlapping.
            /// </summary>
            /// <param name="other">The bin to compare with</param>
            /// <returns>
            /// 0 if the upper bound and lower bound are bit-for-bit equal.
            /// +1 if this bin is lower than the compared bin.
            /// -1 otherwise.
            /// </returns>
            public int CompareTo(Bin other)
            {
                if (UpperBound > other.LowerBound && LowerBound < other.LowerBound)
                {
                    throw new ArgumentException(nameof(other), "The bins cannot be overlapping.");
                }
                if (UpperBound.Equals(other.UpperBound) && LowerBound.Equals(other.LowerBound))
                {
                    return 0;
                }
                if (other.UpperBound <= LowerBound)
                {
                    return 1;
                }
                return -1;
            }

            /// <summary>
            /// Creates a copy of the histogram bin.
            /// </summary>
            /// <returns>A cloned histogram bin.</returns>
            public object Clone()
            {
                return new Bin(LowerBound, UpperBound, Frequency);
            }

            /// <summary>
            /// Checks whether two histogram bins are equal.
            /// </summary>
            /// <returns>True if the bins are equal and false otherwise.</returns>
            public override bool Equals(object obj)
            {
                if (!(obj is Bin))
                {
                    return false;
                }
                Bin bin = (Bin)obj;
                return LowerBound.Equals(bin.LowerBound) && UpperBound.Equals(bin.UpperBound) && Frequency.Equals(bin.Frequency);
            }
        }

        /// <summary>
        /// Constructs a histogram based on the data provided. The Rice Rule is used to set the bin sizes.
        /// </summary>
        /// <param name="data">The data to construct a histogram with.</param>
        public Histogram(IList<double> data)
        {
            // Use the Rice Rule to set the bin count
            NumberOfBins = (int)(Math.Ceiling(2d * Math.Pow(data.Count, 1d / 3d)) + 1d);
            // Get the bin boundaries
            LowerBound = data.Min();
            UpperBound = data.Max();
            BinWidth = (UpperBound - LowerBound) / NumberOfBins;
            // Add bins
            double xl = LowerBound;
            double xu = xl + BinWidth;
            AddBin(new Bin(xl, xu));
            for (int i = 1; i < NumberOfBins; i++)
            {
                xl = xu;
                xu = xl + BinWidth;
                AddBin(new Bin(xl, xu));
            }
            // Guarantee that the upper bound is exactly the data maximum.
            _bins.Last().UpperBound = UpperBound;
            // Add data
            AddData(data);
        }

        /// <summary>
        /// Constructs a histogram with a specific number of bins and the data provided. The histogram limits are derived from the data.
        /// </summary>
        /// <param name="data">The data to construct a histogram with.</param>
        /// <param name="numberOfBins">The number of bins.</param>
        public Histogram(IList<double> data, int numberOfBins)
        {
            // Get the bin boundaries
            NumberOfBins = numberOfBins;
            LowerBound = data.Min();
            UpperBound = data.Max();
            BinWidth = (UpperBound - LowerBound) / NumberOfBins;
            // Add bins
            double xl = LowerBound;
            double xu = xl + BinWidth;
            AddBin(new Bin(xl, xu));
            for (int i = 1; i < NumberOfBins; i++)
            {
                xl = xu;
                xu = xl + BinWidth;
                AddBin(new Bin(xl, xu));
            }
            // Guarantee that the upper bound is exactly the data maximum.
            _bins.Last().UpperBound = UpperBound;
            // Add data
            AddData(data);
        }


        private List<double> _binLimits = new List<double>();
        private List<Bin> _bins = new List<Bin>();
        private bool _areBinsSorted;

        /// <summary>
        /// Get the starting bin value.
        /// </summary>
        public double LowerBound { get; private set; }

        /// <summary>
        /// Get the ending bin value.
        /// </summary>
        public double UpperBound { get; private set; }

        /// <summary>
        /// Returns the number of bins in the histogram.
        /// </summary>
        public int NumberOfBins { get; private set; }

        /// <summary>
        /// Returns the bin width.
        /// </summary>
        public double BinWidth { get; private set; }

        /// <summary>
        /// Returns the length of the data sample.
        /// </summary>
        public int DataCount
        {
            get
            {
                int count = 0;
                for (int i = 0; i < _bins.Count; i++)
                    count += _bins[i].Frequency;
                return count;
            }
        }

        /// <summary>
        /// Get the bin at a specific index.
        /// </summary>
        /// <param name="index">The index of the bin to be returned.</param>
        /// <returns>A copy of the bin.</returns>
        public Bin this[int index]
        {
            get
            {
                SortBins();
                return (Bin)_bins[index].Clone();
            }
        }

        /// <summary>
        /// Gets the mean of the histogram.
        /// </summary>
        public double Mean
        {
            get
            {
                int total = 0;
                double sum = 0d;
                for (int i = 0; i < _bins.Count; i++)
                {
                    sum += _bins[i].Midpoint * _bins[i].Frequency;
                    total += _bins[i].Frequency;
                }
                if (total == 0)
                {
                    return 0d;
                }
                else
                {
                    return sum / total;
                }
            }
        }

        /// <summary>
        /// Gets the median of the histogram.
        /// </summary>
        public double Median
        {
            get
            {
                int total = 0;
                for (int i = 0; i < _bins.Count; i++)
                    total += _bins[i].Frequency;
                int halfTotal = (int)(total / 2d);
                int m = 0;
                int v = 0;
                while (m < _bins.Count)
                {
                    v += _bins[m].Frequency;
                    if (v >= halfTotal)
                    {
                        break;
                    }
                    m += 1;
                }
                return _bins[m].Midpoint;
            }
        }

        /// <summary>
        /// Gets the mode of the histogram.
        /// </summary>
        public double Mode
        {
            get
            {
                int m = 0;
                int curMax = 0;
                for (int i = 0; i < _bins.Count;  i++)
                {
                    if (_bins[i].Frequency > curMax)
                    {
                        curMax = _bins[i].Frequency;
                        m = i;
                    }
                }
                return _bins[m].Midpoint;
            }
        }

        /// <summary>
        /// Gets the standard deviation of the histogram.
        /// </summary>
        public double StandardDeviation
        {
            get
            {
                double stddev = 0d;
                int total = 0;
                double m = Mean;
                for (int i = 0; i < _bins.Count; i++)
                {
                    int vals = _bins[i].Frequency;
                    double diff = _bins[i].Midpoint - m;
                    stddev += diff * diff * vals;
                    total += vals;
                }
                if (total == 0)
                {
                    return 0d;
                }
                else
                {
                    return Math.Sqrt(stddev / total);
                }
            }
        }

        /// <summary>
        /// Add a bin to the bin list.
        /// </summary>
        /// <param name="bin">Histogram bin to add.</param>
        public void AddBin(Bin bin)
        {
            _bins.Add(bin);
            _areBinsSorted = false;
        }

        /// <summary>
        /// Add one data value to the histogram. If the data value falls outside the range of the histogram,
        /// the start or end bin will automatically adapt.
        /// </summary>
        /// <param name="data">The data value to add.</param>
        public void AddData(double data)
        {
            SortBins();
            int index = GetBinIndexOf(data);
            if (data <= LowerBound)
            {
                _bins.First().LowerBound = data;
                _bins.First().Frequency += 1;
                _areBinsSorted = false;
            }
            else if (data >= UpperBound)
            {
                _bins.Last().UpperBound = data;
                _bins.Last().Frequency += 1;
                _areBinsSorted = false;
            }
            else if (index >= 0 && index < NumberOfBins)
            {
                _bins[index].Frequency += 1;
            }
        }

        /// <summary>
        /// Add a sequence of data values to the histogram. If the data value falls outside the range of the histogram,
        /// the start or end bin will automatically adapt.
        /// </summary>
        /// <param name="data">A sequence of data values to add.</param>
        public void AddData(IList<double> data)
        {
            foreach (double d in data)
                AddData(d);
        }

        /// <summary>
        /// Sort the histogram bins if needed.
        /// </summary>
        private void SortBins()
        {
            if (!_areBinsSorted)
            {
                _bins.Sort();
                _binLimits = _bins.Select(x => x.LowerBound).ToList();
                _areBinsSorted = true;
            }
        }

        /// <summary>
        /// Returns the index in the histogram of the bin that contains the value.
        /// </summary>
        /// <param name="value">The value to search for.</param>
        /// <returns>The index of the bin containing the value.</returns>
        public int GetBinIndexOf(double value)
        {
            SortBins();
            if (value < _bins.First().LowerBound || value > _bins.Last().UpperBound)
            {
                throw new ArgumentException("value", "The value is not contained with the histogram limits.");
            }
            int idx = Search.Bisection(value, _binLimits);
            return idx < 0 ? 0 : idx >= _binLimits.Count ? _binLimits.Count - 1 : idx;
        }

    }
}