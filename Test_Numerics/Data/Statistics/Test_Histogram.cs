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
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;

namespace Data.Statistics
{
    /// <summary>
    /// Unit tests for the Histogram class. These methods were tested against values attained from MS Excel or by direct calculation of the 
    /// equation.
    /// </summary>
    /// <remarks>
    /// <para>
    ///      <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     <item>Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_Histogram
    {
        /// <summary>
        /// Test the construction and manipulation of the Bin class used in the Histogram class. 
        /// </summary>
        [TestMethod]
        public void Test_Bin()
        {
            var bin1 = new Histogram.Bin(0, 2);
            var bin2 = new Histogram.Bin(2, 6);
            var bin3 = (Histogram.Bin)bin1.Clone();

            var test1 = bin1.CompareTo(bin2);
            var test2 = bin1.CompareTo(bin3);
            var test3 = bin2.CompareTo(bin3);

            Assert.AreEqual(-1, test1);
            Assert.AreEqual(0, test2);
            Assert.AreEqual(1, test3);

            var test4 = bin1.Equals(bin3);
            var test5 = bin1.Equals(bin2);

            Assert.IsTrue(test4);
            Assert.IsFalse(test5);
        }

        /// <summary>
        /// Test the construction of a basic histogram with no specified number of bins
        /// </summary>
        [TestMethod]
        public void Test_Histogram_Basic()
        {
            var sample1 = new double[] { 122d, 244d, 214d, 173d, 229d, 156d, 212d, 263d, 146d, 183d, 161d, 205d, 135d, 331d, 225d, 174d, 98.8d, 149d, 238d, 262d, 132d, 235d, 216d, 240d, 230d, 192d, 195d, 172d, 173d, 172d, 153d, 142d, 317d, 161d, 201d, 204d, 194d, 164d, 183d, 161d, 167d, 179d, 185d, 117d, 192d, 337d, 125d, 166d, 99.1d, 202d, 230d, 158d, 262d, 154d, 164d, 182d, 164d, 183d, 171d, 250d, 184d, 205d, 237d, 177d, 239d, 187d, 180d, 173d, 174d };
            var hist = new Histogram(sample1);
            int numBins = hist.NumberOfBins;

            // These values come from Microsoft Excel
            var lower = new double[] { 98.8, 122.62, 146.44, 170.26, 194.08, 217.9, 241.72, 265.54, 289.36, 313.18 };
            var upper = new double[] { 122.62, 146.44, 170.26, 194.08, 217.9, 241.72, 265.54, 289.36, 313.18, 337 };
            var frequency = new int[] { 4, 5, 13, 21, 9, 9, 5, 0, 0, 3 };

            for (int i = 0; i < hist.NumberOfBins; i++)
            {
                Assert.AreEqual(hist[i].LowerBound, lower[i], 1E-3);
                Assert.AreEqual(hist[i].UpperBound, upper[i], 1E-3);
                Assert.AreEqual(hist[i].Frequency, frequency[i], 1E-3);
            }
        }

        /// <summary>
        /// Test the construction of a histogram with a specified number of bins
        /// </summary>
        [TestMethod()]
        public void Test_Histogram_Specified()
        {
            var sample1 = new double[] { 122d, 244d, 214d, 173d, 229d, 156d, 212d, 263d, 146d, 183d, 161d, 205d, 135d, 331d, 225d, 174d, 98.8d, 149d, 238d, 262d, 132d, 235d, 216d, 240d, 230d, 192d, 195d, 172d, 173d, 172d, 153d, 142d, 317d, 161d, 201d, 204d, 194d, 164d, 183d, 161d, 167d, 179d, 185d, 117d, 192d, 337d, 125d, 166d, 99.1d, 202d, 230d, 158d, 262d, 154d, 164d, 182d, 164d, 183d, 171d, 250d, 184d, 205d, 237d, 177d, 239d, 187d, 180d, 173d, 174d };

            // These values come from Microsoft Excel
            var lower = new double[] { 98.8d, 110.71d, 122.62d, 134.53d, 146.44d, 158.35d, 170.26d, 182.17d, 194.08d, 205.99d, 217.9d, 229.81d, 241.72d, 253.63d, 265.54d, 277.45d, 289.36d, 301.27d, 313.18d, 325.09d };
            var upper = new double[] { 110.71d, 122.62d, 134.53d, 146.44d, 158.35d, 170.26d, 182.17d, 194.08d, 205.99d, 217.9d, 229.81d, 241.72d, 253.63d, 265.54d, 277.45d, 289.36d, 301.27d, 313.18d, 325.09d, 337.0d };
            var frequency = new int[] { 2, 2, 2, 3, 5, 8, 12, 9, 6, 3, 2, 7, 2, 3, 0, 0, 0, 0, 1, 2 };
            var hist = new Histogram(sample1, 20);

            for (int i = 0; i < hist.NumberOfBins; i++)
            {
                Assert.AreEqual(hist[i].LowerBound, lower[i], 1E-3);
                Assert.AreEqual(hist[i].UpperBound, upper[i], 1E-3);
                Assert.AreEqual(hist[i].Frequency, frequency[i], 1E-3);
            }
        }

        /// <summary>
        /// Test that DataCount functions correctly and returns the number of samples
        /// </summary>
        [TestMethod]
        public void Test_DataCount()
        {
            var sample1 = new double[] { 122d, 244d, 214d, 173d, 229d, 156d, 212d, 263d, 146d, 183d, 161d, 205d, 135d, 331d, 225d, 174d, 98.8d, 149d, 238d, 262d, 132d, 235d, 216d, 240d, 230d, 192d, 195d, 172d, 173d, 172d, 153d, 142d, 317d, 161d, 201d, 204d, 194d, 164d, 183d, 161d, 167d, 179d, 185d, 117d, 192d, 337d, 125d, 166d, 99.1d, 202d, 230d, 158d, 262d, 154d, 164d, 182d, 164d, 183d, 171d, 250d, 184d, 205d, 237d, 177d, 239d, 187d, 180d, 173d, 174d };
            var hist = new Histogram(sample1);

            var test = hist.DataCount;
            Assert.AreEqual(sample1.Count(), test);
        }

        /// <summary>
        /// Test the statistics of the histogram: mean, median, mode, and standard deviation.
        /// </summary>
        [TestMethod]
        public void Test_Stats()
        {
            var sample1 = new double[] { 122d, 244d, 214d, 173d, 229d, 156d, 212d, 263d, 146d, 183d, 161d, 205d, 135d, 331d, 225d, 174d, 98.8d, 149d, 238d, 262d, 132d, 235d, 216d, 240d, 230d, 192d, 195d, 172d, 173d, 172d, 153d, 142d, 317d, 161d, 201d, 204d, 194d, 164d, 183d, 161d, 167d, 179d, 185d, 117d, 192d, 337d, 125d, 166d, 99.1d, 202d, 230d, 158d, 262d, 154d, 164d, 182d, 164d, 183d, 171d, 250d, 184d, 205d, 237d, 177d, 239d, 187d, 180d, 173d, 174d };
            var hist = new Histogram(sample1);
            var lower = new double[] { 98.8, 122.62, 146.44, 170.26, 194.08, 217.9, 241.72, 265.54, 289.36, 313.18 };
            var upper = new double[] { 122.62, 146.44, 170.26, 194.08, 217.9, 241.72, 265.54, 289.36, 313.18, 337 };
            var frequency = new int[] { 4, 5, 13, 21, 9, 9, 5, 0, 0, 3 };
            
            var products = new double[hist.NumberOfBins];
            // basic formula for mean
            for(int i = 0; i < hist.NumberOfBins; i++)
            {
                double mid = (upper[i] + lower[i]) / 2;
                products[i] = mid*frequency[i];
            }
            var trueMean = products.Sum() / frequency.Sum();
            double mean = hist.Mean;

            Assert.AreEqual(trueMean, mean, 1E-6);


            int total = frequency.Sum();
            int halfTotal = total / 2;
            int m = 0; int n = 0;
            // basic formula for median
            while (m < hist.NumberOfBins)
            {
                n += frequency[m];
                if (n >= halfTotal) { break; }
                m++;
            }
            double trueMedian = (upper[m] + lower[m]) / 2; 
            double median = hist.Median;

            Assert.AreEqual(trueMedian, median, 1E-6);


            double max = frequency.Max();
            int index = 0;
            // basic formula for mode
            for(int i = 0; i < frequency.Length; i++)
            {
                if(frequency[i] == max) {index = i;  break;}
            }
            double trueMode = (upper[index] + lower[index]) / 2;
            double mode = hist.Mode;

            Assert.AreEqual(trueMode, mode, 1E-6);


            double dev = 0; int tot = 0; double mu = trueMean; 
            // basic formula for standard deviation
            for(int i = 0; i < hist.NumberOfBins; i++)
            {
                dev += (((upper[i] + lower[i]) / 2) - mu) * (((upper[i] + lower[i]) / 2) - mu) * frequency[i];
                tot += frequency[i];
            }
            double trueSD = Math.Sqrt(dev / tot);
            double sd = hist.StandardDeviation;

            Assert.AreEqual(trueSD, sd, 1E-6);
        }

        /// <summary>
        /// Test the indexing method of the histogram
        /// </summary>
        [TestMethod]
        public void Test_Indexing()
        {
            var sample1 = new double[] { 122d, 244d, 214d, 173d, 229d, 156d, 212d, 263d, 146d, 183d, 161d, 205d, 135d, 331d, 225d, 174d, 98.8d, 149d, 238d, 262d, 132d, 235d, 216d, 240d, 230d, 192d, 195d, 172d, 173d, 172d, 153d, 142d, 317d, 161d, 201d, 204d, 194d, 164d, 183d, 161d, 167d, 179d, 185d, 117d, 192d, 337d, 125d, 166d, 99.1d, 202d, 230d, 158d, 262d, 154d, 164d, 182d, 164d, 183d, 171d, 250d, 184d, 205d, 237d, 177d, 239d, 187d, 180d, 173d, 174d };
            var hist = new Histogram(sample1);
            for(int i = 0; i < hist.NumberOfBins; i++)
            {
                int index = hist.GetBinIndexOf(hist[i].Midpoint);
                Assert.AreEqual(i, index);
            }

            var hist2 = new Histogram(sample1, 20);
            for (int i = 0; i < hist2.NumberOfBins; i++)
            {
                int index = hist2.GetBinIndexOf(hist2[i].Midpoint);
                Assert.AreEqual(i, index);
            }
        }
    }
}
