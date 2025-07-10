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

namespace Numerics.Data
{
    /// <summary>
    /// A base class for interpolation.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// Interpolation is an estimation method to find new data points within the range of a discrete set of known data points.
    /// </para>
    /// <b> References: </b>
    /// <see href="https://en.wikipedia.org/wiki/Interpolation"/>
    /// </remarks>
    [Serializable]
    public abstract class Interpolater
    {

        /// <summary>
        /// Construct a new interpolation class. 
        /// </summary>
        /// <param name="xValues">List of x-values. The x-Values must be monotonic, either increasing or decreasing.</param>
        /// <param name="yValues">List of y-values.</param>
        /// <param name="sortOrder">The sort order of the x-values, either ascending or descending. Default = Ascending. </param>
        public Interpolater(IList<double> xValues, IList<double> yValues, SortOrder sortOrder = SortOrder.Ascending)
        {
            Count = xValues.Count;

            // Validate 
            if (yValues.Count != Count) throw new ArgumentException(nameof(xValues), "The x and y lists must be the same length.");
            if (Count < 2) throw new ArgumentException(nameof(xValues), "The x list is too small. It must have at least 2 values.");
            for (int i = 1; i < xValues.Count; ++i)
            {
                if (xValues[i] == xValues[i - 1]) throw new ArgumentException(nameof(xValues), "All x values should be unique.");
                if (sortOrder == SortOrder.Ascending && xValues[i] < xValues[i - 1]) throw new ArgumentException(nameof(xValues), "The x values are not in ascending order.");
                if (sortOrder == SortOrder.Descending && xValues[i] > xValues[i - 1]) throw new ArgumentException(nameof(xValues), "The x values are not in descending order.");
            }
            this.XValues = xValues;
            this.YValues = yValues;     
            deltaStart = Math.Min(1, (int)Math.Pow((double)Count, 0.25));
            SortOrder = sortOrder;
            
        }

        /// <summary>
        /// The size of the xy-value list.
        /// </summary>
        public int Count { get; private set; } = 0;

        /// <summary>
        /// Keeps track of the starting search location. 
        /// </summary>
        public int SearchStart { get; set; } = 0;

        /// <summary>
        /// Keeps track of the difference is start locations. 
        /// </summary>
        protected int deltaStart = 0;

        /// <summary>
        /// Determines which search method to use. If values are correlated, use the Hunt method. 
        /// </summary>
        protected bool correlated = false;

        /// <summary>
        /// Determines whether to use a smart searching algorithm or just sequential search.
        /// </summary>
        public bool UseSmartSearch { get; set; } = true;

        /// <summary>
        /// The array of x-values.
        /// </summary>
        public IList<double> XValues { get; protected set; }

        /// <summary>
        /// The array of y-values.
        /// </summary>
        public IList<double> YValues { get; protected set; }

        /// <summary>
        /// The sort order of the x-values, either ascending or descending. Default = Ascending. 
        /// </summary>
        public SortOrder SortOrder { get; protected set; }

        /// <summary>
        /// Given a value x, returns an interpolated value.
        /// </summary>
        /// <param name="x">The value to interpolate.</param>
        /// <param name="index">The zero-based index of the lower bound of the interpolation interval.</param>
        public abstract double BaseInterpolate(double x, int index);

        /// <summary>
        /// Given a list of x-values, returns an array of interpolated values.
        /// </summary>
        /// <param name="x">The list of values to interpolate.</param>
        public double[] Interpolate(IList<double> x)
        {
            var values = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                values[i] = Interpolate(x[i]);
            return values;
        }

        /// <summary>
        /// Given a value x, returns an interpolated value. 
        /// </summary>
        /// <param name="x">The value to interpolate.</param>
        public double Interpolate(double x)
        {        
            return BaseInterpolate(x, Search(x));
        }

        /// <summary>
        /// Search for the lower bound of the interpolation interval. This method updates whether the values being searched on repeated calls are correlated, 
        /// and saves search value for future use on the next call.  
        /// </summary>
        /// <param name="x">The value to search for.</param>
        /// <returns>
        /// The lower bound of the interpolation interval
        /// </returns>
        public int Search(double x)
        {
            int start = 0;
            if (UseSmartSearch)
            {
                start = correlated ? HuntSearch(x) : BisectionSearch(x);
            }
            else
            {
                start = SequentialSearch(x);
            }
            correlated = Math.Abs(start - SearchStart) > deltaStart ? false : true;
            SearchStart = start < 0 || start >= Count ? 0 : start;
            return start;
        }

        /// <summary>
        /// Searches for the lower bound of the location of a value using a sequential search method.
        /// </summary>
        /// <param name="x">The value to search for.</param>
        /// <returns>
        /// The lower bound of the location of a value
        /// </returns>
        public int SequentialSearch(double x)
        {
            int jl = SearchStart;
            if ((SortOrder == SortOrder.Ascending && x < XValues[0]) ||
                (SortOrder == SortOrder.Descending && x > XValues[0]))
            {
                return 0;
            }
            else if ((SortOrder == SortOrder.Ascending && x > XValues[Count - 1]) ||
                        (SortOrder == SortOrder.Descending && x < XValues[Count - 1]))
            {
                return Count - 2;
            }
            else if ((SortOrder == SortOrder.Ascending && x < XValues[SearchStart]) ||
                        (SortOrder == SortOrder.Descending && x > XValues[SearchStart]))
            {
                jl = 0;
            }
            for (int i = jl; i < Count; i++)
            {
                if ((SortOrder == SortOrder.Ascending && x <= XValues[i]) ||
                    (SortOrder == SortOrder.Descending && x >= XValues[i]))
                {               
                    jl = i - 1;
                    break;
                }
            }
            return jl;
        }

        /// <summary>
        /// Searches for the lower bound of the location of a value using a bisection search method.
        /// </summary>
        /// <param name="x">The value to search for.</param>
        /// <returns>
        /// The lower bound of the location of a value
        /// </returns>
        public int BisectionSearch(double x)
        {
            int ju = Count - 1, jm, jl = 0;
            bool ascnd = SortOrder == SortOrder.Ascending;
            while (ju - jl > 1)
            {
                jm = (ju + jl) >> 1;
                if (x >= XValues[jm] == ascnd)
                    jl = jm;
                else
                    ju = jm;
            }
            return jl;
        }

        /// <summary>
        /// Searches for the lower bound of the location of a value using a hunt search method.
        /// </summary>
        /// <param name="x">The value to search for.</param>
        /// <returns>
        /// The lower bound of the location of a value
        /// </returns>
        public int HuntSearch(double x)
        {
            int jl = SearchStart, jm, ju, inc = 1;
            bool ascnd = SortOrder == SortOrder.Ascending;
            if (jl < 0 || jl > Count - 1)
            {
                jl = 0;
                ju = Count - 1;
            }
            else
            {
                if (x >= XValues[jl] == ascnd)
                {
                    for (; ; )
                    {
                        ju = jl + inc;
                        if (ju >= Count - 1) { ju = Count - 1; break; }
                        else if (x < XValues[ju] == ascnd) break;
                        else
                        {
                            jl = ju;
                            inc += inc;
                        }
                    }
                }
                else
                {
                    ju = jl;
                    for (; ; )
                    {
                        jl = jl - inc;
                        if (jl <= 0) { jl = 0; break; }
                        else if (x >= XValues[jl] == ascnd) break;
                        else
                        {
                            ju = jl;
                            inc += inc;
                        }
                    }
                }
            }
            while (ju - jl > 1)
            {
                jm = (ju + jl) >> 1;
                if (x >= XValues[jm] == ascnd)
                    jl = jm;
                else
                    ju = jm;
            }
            return jl;
        }

    }
}
