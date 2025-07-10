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
    /// A class for searching a list.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// This class has 3 search methods: Sequential, Bisection, and Hunt.
    /// <list type="bullet">
    /// <item>
    /// <description>
    /// Sequential search, or linear search, checks each element in a list until a match is found or the whole list has been searched.
    /// </description>
    /// </item>
    /// <item>
    /// <description> 
    /// The bisection method, also known as binary search, starts off by comparing the target value (i.e. value being searched for in the array) 
    /// against the middle element of the list. If the target value and middle element are not equal, the half of the list in which the target 
    /// cannot lie is eliminated. If the list  sorted in ascending order, if the target is less than the middle element, the second half of the list of values 
    /// greater than the middle element is eliminated. If the target is greater than the middle element then the first half of the list is eliminated. Then, the algorithm starts
    /// again, comparing the middle element of the new, smaller list against the target element. This continues until the target value is found or until the remaining half 
    /// being empty.
    /// </description>
    /// </item>
    /// </list>
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/Linear_search"/>
    /// </description></item>
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/Binary_search"/>
    /// </description></item>
    /// </list>
    /// </remarks>
    [Serializable]
    public class Search
    {
        /// <summary>
        /// Enumeration of search methods.
        /// </summary>
        public enum Method
        {
            /// <summary>
            /// Search sequentially.
            /// </summary>
            Sequential,

            /// <summary>
            /// Search with bisection method.
            /// </summary>
            Bisection,

            /// <summary>
            /// Search using the hunt algorithm.
            /// </summary>
            Hunt
        }

        #region Sequential

        /// <summary>
        /// Searches for the lower bound of the location of a value using a sequential search method.
        /// </summary>
        /// <param name="x">The value to search for.</param>
        /// <param name="values">Array of values to search within.</param>
        /// <param name="start">Optional. Location to start the search of the arrays. Default = 0.</param>
        /// <param name="order">Optional. The sort order of the values, either ascending or descending. Default = Ascending.</param>
        /// <returns>
        /// The lower bound of the location to be used in interpolation. E.g., X1 is the lower, X2 is the upper (X1 + 1).
        /// </returns>
        public static int Sequential(double x, IList<double> values, int start = 0, SortOrder order = SortOrder.Ascending)
        {
            // variables
            int N = values.Count;

            // validate
            if (start < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(start), "The search starting point must be non-negative.");
            }
            else if (start >= N)
            {
                throw new ArgumentOutOfRangeException(nameof(start), "The search starting point cannot be greater than the length of the X array.");
            }

            // check if the order is ascending, and if the x value is outside of the range of the values
            if (order == SortOrder.Ascending)
            {
                if (x < values[0])
                {
                    return -1;
                }
                else if (x == values[0])
                {
                    return 0;
                }
                else if (x == values[N - 1])
                {
                    return N - 1;
                }
                else if (x > values[N - 1])
                {
                    return N;
                }

                for (int i = start; i < N; i++)
                {
                    if (x <= values[i])
                    {
                        return i - 1;
                    }
                }

            }
            else
            {

                if (x > values[0])
                {
                    return -1;
                }
                else if (x == values[0])
                {
                    return 0;
                }
                else if (x == values[N - 1])
                {
                    return N - 1;
                }
                else if (x < values[N - 1])
                {
                    return N;
                }

                for (int i = start; i < N; i++)
                {
                    if (x >= values[i])
                    {
                        return i - 1;
                    }
                }
            }

            return 0;
        }

        /// <summary>
        /// Searches for the lower bound of the location of a value using a sequential search method.
        /// </summary>
        /// <param name="x">The value to search for.</param>
        /// <param name="orderedPairedData">Ordered paired data to search within.</param>
        /// <param name="start">Optional. Location to start the search of the arrays. Default = 0.</param>
        /// <returns>
        /// The lower bound of the location to be used in interpolation. E.g., X1 is the lower, X2 is the upper (X1 + 1).
        /// </returns>
        public static int Sequential(double x, OrderedPairedData orderedPairedData, int start = 0)
        {
            // variables
            int N = orderedPairedData.Count;

            // validate
            if (start < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(start), "The search starting point must be non-negative.");
            }
            else if (start >= N)
            {
                throw new ArgumentOutOfRangeException(nameof(start), "The search starting point cannot be greater than the length of the X array.");
            }

            // check if the order is ascending, and if the x value is outside of the range of the values
            if (orderedPairedData.OrderX == SortOrder.Ascending)
            {
                if (x < orderedPairedData[0].X)
                {
                    return -1;
                }
                else if (x == orderedPairedData[0].X)
                {
                    return 0;
                }
                else if (x == orderedPairedData[N - 1].X)
                {
                    return N - 1;
                }
                else if (x > orderedPairedData[N - 1].X)
                {
                    return N;
                }

                for (int i = start; i < N; i++)
                {
                    if (x <= orderedPairedData[i].X)
                    {
                        return i - 1;
                    }
                }

            }
            else
            {

                if (x > orderedPairedData[0].X)
                {
                    return -1;
                }
                else if (x == orderedPairedData[0].X)
                {
                    return 0;
                }
                else if (x == orderedPairedData[N - 1].X)
                {
                    return N - 1;
                }
                else if (x < orderedPairedData[N - 1].X)
                {
                    return N;
                }

                for (int i = start; i < N; i++)
                {
                    if (x >= orderedPairedData[i].X)
                    {
                        return i - 1;
                    }
                }
            }

            return 0;

        }

        /// <summary>
        /// Searches for the lower bound of the location of a value using a sequential search method.
        /// </summary>
        /// <param name="x">The value to search for.</param>
        /// <param name="ordinates">List of ordinates to search within.</param>
        /// <param name="start">Optional. Location to start the search of the arrays. Default = 0.</param>
        /// <param name="order">Optional. The sort order of the values, either ascending or descending. Default = Ascending.</param>
        /// <returns>
        /// The lower bound of the location to be used in interpolation. E.g., X1 is the lower, X2 is the upper (X1 + 1).
        /// </returns>
        public static int Sequential(double x, IList<Ordinate> ordinates, int start = 0, SortOrder order = SortOrder.Ascending)
        {
            // variables
            int N = ordinates.Count;

            // validate
            if (start < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(start), "The search starting point must be non-negative.");
            }
            else if (start >= N)
            {
                throw new ArgumentOutOfRangeException(nameof(start), "The search starting point cannot be greater than the length of the X array.");
            }

            // check if the order is ascending, and if the x value is outside of the range of the values
            if (order == SortOrder.Ascending)
            {
                if (x < ordinates[0].X)
                {
                    return -1;
                }
                else if (x == ordinates[0].X)
                {
                    return 0;
                }
                else if (x == ordinates[N - 1].X)
                {
                    return N - 1;
                }
                else if (x > ordinates[N - 1].X)
                {
                    return N;
                }

                for (int i = start; i < N; i++)
                {
                    if (x <= ordinates[i].X)
                    {
                        return i - 1;
                    }
                }

            }
            else
            {

                if (x > ordinates[0].X)
                {
                    return -1;
                }
                else if (x == ordinates[0].X)
                {
                    return 0;
                }
                else if (x == ordinates[N - 1].X)
                {
                    return N - 1;
                }
                else if (x < ordinates[N - 1].X)
                {
                    return N;
                }

                for (int i = start; i < N; i++)
                {
                    if (x >= ordinates[i].X)
                    {
                        return i - 1;
                    }
                }
            }

            return 0;
        }

        #endregion

        #region Bisection

        /// <summary>
        /// Searches for the lower bound of the location of a value using a bisection search method.
        /// </summary>
        /// <param name="x">The value to search for.</param>
        /// <param name="values">Array of values to search within.</param>
        /// <param name="start">Optional. Location to start the search of the arrays. Default = 0.</param>
        /// <param name="order">Optional. The sort order of the values, either ascending or descending. Default = Ascending.</param>
        /// <returns>
        /// The lower bound of the location to be used in interpolation. E.g., X1 is the lower, X2 is the upper (X1 + 1).
        /// </returns>
        public static int Bisection(double x, IList<double> values, int start = 0, SortOrder order = SortOrder.Ascending)
        {
            // variables
            int N = values.Count;

            // validate
            if (start < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(start), "The search starting point must be non-negative.");
            }
            else if (start >= N)
            {
                throw new ArgumentOutOfRangeException(nameof(start), "The search starting point cannot be greater than the length of the X array.");
            }

            // check if the order is ascending, if the X is outside of the range of the values
            if (order == SortOrder.Ascending)
            {
                if (x < values[0])
                {
                    return -1;
                }
                else if (x == values[0])
                {
                    return 0;
                }
                else if (x == values[N - 1])
                {
                    return N - 1;
                }
                else if (x > values[N - 1])
                {
                    return N;
                }
            }
            else
            {
                if (x > values[0])
                {
                    return -1;
                }
                else if (x == values[0])
                {
                    return 0;
                }
                else if (x == values[N - 1])
                {
                    return N - 1;
                }
                else if (x < values[N - 1])
                {
                    return N;
                }
            }

            // Perform bisection search
            int xlo = start, xhi = N, xm = 0;
            while (xhi - xlo > 1)
            {
                xm = xlo + (xhi - xlo >> 1); 
                if (x >= values[xm] && order == SortOrder.Ascending)
                {
                    xlo = xm;
                }
                else
                {
                    xhi = xm;
                }
            }
            // Return XLO
            return xlo;
        }

        /// <summary>
        /// Searches for the lower bound of the location of a value using a bisection search method.
        /// </summary>
        /// <param name="x">The value to search for.</param>
        /// <param name="orderedPairedData">Ordered paired data to search within.</param>
        /// <param name="start">Optional. Location to start the search of the arrays. Default = 0.</param>
        /// <returns>
        /// The lower bound of the location to be used in interpolation. E.g., X1 is the lower, X2 is the upper (X1 + 1).
        /// </returns>
        public static int Bisection(double x, OrderedPairedData orderedPairedData, int start = 0)
        {
            // variables
            int N = orderedPairedData.Count;

            // validate
            if (start < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(start), "The search starting point must be non-negative.");
            }
            else if (start >= N)
            {
                throw new ArgumentOutOfRangeException(nameof(start), "The search starting point cannot be greater than the length of the X array.");
            }

            // check if the order is ascending, if the X is outside of the range of the array
            if (orderedPairedData.OrderX == SortOrder.Ascending)
            {
                if (x < orderedPairedData[0].X)
                {
                    return -1;
                }
                else if (x == orderedPairedData[0].X)
                {
                    return 0;
                }
                else if (x == orderedPairedData[N - 1].X)
                {
                    return N - 1;
                }
                else if (x > orderedPairedData[N - 1].X)
                {
                    return N;
                }
            }
            else
            {
                if (x > orderedPairedData[0].X)
                {
                    return -1;
                }
                else if (x == orderedPairedData[0].X)
                {
                    return 0;
                }
                else if (x == orderedPairedData[N - 1].X)
                {
                    return N - 1;
                }
                else if (x < orderedPairedData[N - 1].X)
                {
                    return N;
                }
            }

            // Perform bisection search
            int xlo = start, xhi = N, xm = 0;
            while (xhi - xlo > 1)
            {
                xm = xlo + (xhi - xlo >> 1);
                if (x >= orderedPairedData[xm].X && orderedPairedData.OrderX == SortOrder.Ascending)
                {
                    xlo = xm;
                }
                else
                {
                    xhi = xm;
                }
            }
            // Return XLO
            return xlo;
        }

        /// <summary>
        /// Searches for the lower bound of the location of a value using a bisection search method.
        /// </summary>
        /// <param name="x">The value to search for.</param>
        /// <param name="ordinates">List of ordinates to search within.</param>
        /// <param name="start">Optional. Location to start the search of the arrays. Default = 0.</param>
        /// <param name="order">Optional. The sort order of the values, either ascending or descending. Default = Ascending.</param>
        /// <returns>
        /// The lower bound of the location to be used in interpolation. E.g., X1 is the lower, X2 is the upper (X1 + 1).
        /// </returns>
        public static int Bisection(double x, IList<Ordinate> ordinates, int start = 0, SortOrder order = SortOrder.Ascending)
        {
            // variables
            int N = ordinates.Count;

            // validate
            if (start < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(start), "The search starting point must be non-negative.");
            }
            else if (start >= N)
            {
                throw new ArgumentOutOfRangeException(nameof(start), "The search starting point cannot be greater than the length of the X array.");
            }

            // check if the order is ascending, if the X is outside of the range of the array
            if (order == SortOrder.Ascending)
            {
                if (x < ordinates[0].X)
                {
                    return -1;
                }
                else if (x == ordinates[0].X)
                {
                    return 0;
                }
                else if (x == ordinates[N - 1].X)
                {
                    return N - 1;
                }
                else if (x > ordinates[N - 1].X)
                {
                    return N;
                }
            }
            else
            {
                if (x > ordinates[0].X)
                {
                    return -1;
                }
                else if (x == ordinates[0].X)
                {
                    return 0;
                }
                else if (x == ordinates[N - 1].X)
                {
                    return N - 1;
                }
                else if (x < ordinates[N - 1].X)
                {
                    return N;
                }
            }

            // Perform bisection search
            int xlo = start, xhi = N, xm = 0;
            while (xhi - xlo > 1)
            {
                xm = xlo + (xhi - xlo >> 1);
                if (x >= ordinates[xm].X && order == SortOrder.Ascending)
                {
                    xlo = xm;
                }
                else
                {
                    xhi = xm;
                }
            }
            // Return XLO
            return xlo;
        }

        #endregion

        #region Hunt

        /// <summary>
        /// Searches for the lower bound of the location of a value using the Hunt method.
        /// </summary>
        /// <param name="x">The value to search for.</param>
        /// <param name="values">Array of values.</param>
        /// <param name="start">Optional. Location to start the search of the arrays. Default = 0.</param>
        /// <param name="order">Optional. The sort order of the values, either ascending or descending. Default = Ascending.</param>
        /// <returns>
        /// The lower bound of the location to be used in interpolation. E.g., X1 is the lower, X2 is the upper (X1 + 1).
        /// </returns>
        public static int Hunt(double x, IList<double> values, int start = 0, SortOrder order = SortOrder.Ascending)
        {
            int N = values.Count;
            
            // validate
            if (start < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(start), "The search starting point must be non-negative.");
            }
            else if (start >= N)
            {
                throw new ArgumentOutOfRangeException(nameof(start), "The search starting point cannot be greater than the length of the X array.");
            }

            // check if the order is ascending, if the X is outside of the range of the values
            if (order == SortOrder.Ascending)
            {
                if (x < values[0])
                {
                    return -1;
                }
                else if (x == values[0])
                {
                    return 0;
                }
                else if (x == values[N - 1])
                {
                    return N - 1;
                }
                else if (x > values[N - 1])
                {
                    return N;
                }
            }
            else
            {
                if (x > values[0])
                {
                    return -1;
                }
                else if (x == values[0])
                {
                    return 0;
                }
                else if (x == values[N - 1])
                {
                    return N - 1;
                }
                else if (x < values[N - 1])
                {
                    return N;
                }
            }

            // XLO is defined by the search start
            int xlo = start, xhi = N, xm = 0, inc = 1;

            // Perform the hunt search algorithm
            if (xlo < 0 || xlo > N - 1)
            {
                // The input guess is not useful.
                // Go immediately to bisection. 
                xlo = 0;
                xhi = N + 1;
            }
            else
            {
                if (x >= values[xlo] && order == SortOrder.Ascending)
                {
                    // Hunt up
                    for (;;)
                    {
                        // Not done hunting so double the increment
                        xhi = xlo + inc;
                        if (xhi >= N - 1) { xhi = N - 1; break; }
                        else if (x < values[xhi] && order == SortOrder.Ascending) break;
                        else
                        {
                            xlo = xhi;
                            inc += inc;
                        }
                    }
                }
                else
                {
                    // Hunt down
                    xhi = xlo;               
                    for (;;)
                    {
                        xlo = xlo - inc;
                        if (xlo <= 0) { xlo = 0; break; }
                        else if (x >= values[xlo] && order == SortOrder.Ascending) break;
                        else
                        {
                            xhi = xlo;
                            inc += inc;
                        }
                    }
                }
            }

            // Hunt is done, so begin the final bisection
            while (xhi - xlo > 1)
            {
                xm = xlo + (xhi - xlo >> 1);
                if (x >= values[xm] && order == SortOrder.Ascending)
                {
                    xlo = xm;
                }
                else
                {
                    xhi = xm;
                }
            }
            // Return XLO
            return xlo;
        }

        /// <summary>
        /// Searches for the lower bound of the location of a value using the HUNT method.
        /// </summary>
        /// <param name="xValue">The value to search for.</param>
        /// <param name="orderedPairedData">Ordered paired data.</param>
        /// <param name="start">Optional. Location to start the search of the arrays. Default = 0.</param>
        /// <returns>
        /// The lower bound of the location to be used in interpolation. E.g., X1 is the lower, X2 is the upper (X1 + 1).
        /// </returns>
        /// <remarks>
        /// Given an array of X and given a value of X, returns a location J such that X is between XArray(J) and XArray(J+1). The XArray must be
        /// monotonic, either increasing or decreasing. J = 0 or J = XArray.Length is returned to indicate that X is out of range. The SearchStart is
        /// taken as the initial guess for J on output.
        /// </remarks>
        public static int Hunt(double xValue, OrderedPairedData orderedPairedData, int start = 0)
        {
            int N = orderedPairedData.Count - 1;
            double X = xValue;
            int XLO = start;
            int XHI;
            int XM;
            var ASCND = default(bool);
            int INC;

            // validate
            if (start < 0)
            {
                throw new ArgumentOutOfRangeException("searchStart", "The search starting point must be non-negative.");
            }
            else if (start > N)
            {
                throw new ArgumentOutOfRangeException("searchStart", "The search starting point cannot be greater than the length of the X array.");
            }

            // check if the order is ascending, if the X is outside of the range of the array
            if (orderedPairedData.OrderX == SortOrder.Ascending)
            {
                if (xValue < orderedPairedData[0].X)
                {
                    return -1;
                }
                else if (xValue == orderedPairedData[0].X)
                {
                    return 0;
                }
                else if (xValue == orderedPairedData[N].X)
                {
                    return N;
                }
                else if (xValue > orderedPairedData[N].X)
                {
                    return N + 1;
                }
            }
            else if (xValue > orderedPairedData[0].X)
            {
                return -1;
            }
            else if (xValue == orderedPairedData[0].X)
            {
                return 0;
            }
            else if (xValue == orderedPairedData[N].X)
            {
                return N;
            }
            else if (xValue < orderedPairedData[N].X)
            {
                return N + 1;
            }

            // XLO is defined by the search start
            XLO = start;

            // Perform the hunt search algorithm
            if (XLO <= 0 | XLO > N)
            {
                // The input guess is not useful.
                // Go immediately to bisection. 
                XLO = 0;
                XHI = N + 1;
            }
            else
            {
                // Set the hunting increment. 
                INC = 1;
                if (X >= orderedPairedData[XLO].X == ASCND)
                {
                    // Hunt up
                    XHI = XLO + 1;
                    while (X >= orderedPairedData[XHI].X == ASCND)
                    {
                        // Not done hunting so double the increment
                        XLO = XHI;
                        INC += INC;
                        XHI = XLO + INC;
                        if (XHI > N)
                        {
                            // Done hunting since of end of table
                            XHI = N + 1;
                            break;
                        }
                    }
                }
                // Done hunting, value bracketed. 
                else
                {
                    // Hunt down
                    XHI = XLO;
                    XHI = XLO - 1;
                    while (X < orderedPairedData[XLO].X == ASCND)
                    {
                        // Not done hunting so double the increment
                        XHI = XLO;
                        INC += INC;
                        XLO = XHI - INC;
                        if (XLO < 1)
                        {
                            // Done hunting since of end of table
                            XLO = 0;
                            break;
                        }
                    }
                    // Done hunting, value bracketed. 
                }
            }
            // Hunt is done, so begin the final bisection
            while (XHI - XLO > 1)
            {
                XM = XLO + (XHI - XLO >> 1);
                if (X >= orderedPairedData[XM].X && ASCND == true)
                {
                    XLO = XM;
                }
                else
                {
                    XHI = XM;
                }
            }
            // Return XLO
            return XLO;
        }

        /// <summary>
        /// Searches for the lower bound of the location of a value using the HUNT method.
        /// </summary>
        /// <param name="xValue">The value to search for.</param>
        /// <param name="ordinateData">The list of ordinates.</param>
        /// <param name="start">Optional. Location to start the search of the arrays. Default = 0.</param>
        /// <param name="order">Optional. Ascending or descending. Default = ASC.</param>
        /// <returns>
        /// The lower bound of the location to be used in interpolation. E.g., X1 is the lower, X2 is the upper (X1 + 1).
        /// </returns>
        /// <remarks>
        /// Given an array of X and given a value of X, returns a location J such that X is between XArray(J) and XArray(J+1). The XArray must be
        /// monotonic, either increasing or decreasing. J = 0 or J = XArray.Length is returned to indicate that X is out of range. The SearchStart is
        /// taken as the initial guess for J on output.
        /// </remarks>
        public static int Hunt(double xValue, IList<Ordinate> ordinateData, int start = 0, SortOrder order = SortOrder.Ascending)
        {
            int N = ordinateData.Count - 1;
            double X = xValue;
            int XLO = start;
            int XHI;
            int XM;
            bool ASCND;
            int INC;

            // validate
            if (start < 0)
            {
                throw new ArgumentOutOfRangeException("searchStart", "The search starting point must be non-negative.");
            }
            else if (start > N)
            {
                throw new ArgumentOutOfRangeException("searchStart", "The search starting point cannot be greater than the length of the X array.");
            }

            // check if the order is ascending, if the X is outside of the range of the array
            if (order == SortOrder.Ascending)
            {
                ASCND = true;
                if (xValue < ordinateData[0].X)
                {
                    return -1;
                }
                else if (xValue == ordinateData[0].X)
                {
                    return 0;
                }
                else if (xValue == ordinateData[N].X)
                {
                    return N;
                }
                else if (xValue > ordinateData[N].X)
                {
                    return N + 1;
                }
            }
            else
            {
                ASCND = false;
                if (xValue > ordinateData[0].X)
                {
                    return -1;
                }
                else if (xValue == ordinateData[0].X)
                {
                    return 0;
                }
                else if (xValue == ordinateData[N].X)
                {
                    return N;
                }
                else if (xValue < ordinateData[N].X)
                {
                    return N + 1;
                }
            }

            // XLO is defined by the search start
            XLO = start;

            // Perform the hunt search algorithm
            if (XLO <= 0 | XLO > N)
            {
                // The input guess is not useful.
                // Go immediately to bisection. 
                XLO = 0;
                XHI = N + 1;
            }
            else
            {
                // Set the hunting increment. 
                INC = 1;
                if (X >= ordinateData[XLO].X == ASCND)
                {
                    // Hunt up
                    XHI = XLO + 1;
                    while (X >= ordinateData[XHI].X == ASCND)
                    {
                        // Not done hunting so double the increment
                        XLO = XHI;
                        INC += INC;
                        XHI = XLO + INC;
                        if (XHI > N)
                        {
                            // Done hunting since of end of table
                            XHI = N + 1;
                            break;
                        }
                    }
                }
                // Done hunting, value bracketed. 
                else
                {
                    // Hunt down
                    XHI = XLO;
                    XHI = XLO - 1;
                    while (X < ordinateData[XLO].X == ASCND)
                    {
                        // Not done hunting so double the increment
                        XHI = XLO;
                        INC += INC;
                        XLO = XHI - INC;
                        if (XLO < 1)
                        {
                            // Done hunting since of end of table
                            XLO = 0;
                            break;
                        }
                    }
                    // Done hunting, value bracketed. 
                }
            }
            // Hunt is done, so begin the final bisection
            while (XHI - XLO > 1)
            {
                XM = XLO + (XHI - XLO >> 1);
                if (X >= ordinateData[XM].X && ASCND == true)
                {
                    XLO = XM;
                }
                else
                {
                    XHI = XM;
                }
            }
            // Return XLO
            return XLO;
        }

        #endregion

    }
}
