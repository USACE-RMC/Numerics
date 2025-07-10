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
using Numerics.Distributions;

namespace Numerics.Data
{
    /// <summary>
    /// Class to store ordinate information where X and Y are stored as double precision numbers.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    /// <list type="bullet">
    /// <item><description>
    ///     Woodrow Fields, USACE Risk Management Center, woodrow.l.fields@usace.army.mil 
    /// </description></item>
    /// <item><description>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil 
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public struct Ordinate
    {

        #region Construction

        /// <summary>
        /// X-Y ordinate.
        /// </summary>
        /// <param name="xValue">The x-value.</param>
        /// <param name="yValue">The y-value.</param>
        public Ordinate(double xValue, double yValue)
        {
            X = xValue;
            Y = yValue;
            IsValid = true;
            if (double.IsInfinity(X) | double.IsNaN(X) || double.IsInfinity(Y) | double.IsNaN(Y)) { IsValid = false; }
        }


        /// <summary>
        /// Constructs new ordinate from XElement.
        /// </summary>
        /// <param name="xElement">The XElement to deserialize.</param>
        public Ordinate(XElement xElement)
        {
            double x = 0, y = 0;
            if (xElement.Attribute(nameof(X)) != null) double.TryParse(xElement.Attribute(nameof(X)).Value, NumberStyles.Any, CultureInfo.InvariantCulture, out x);
            if (xElement.Attribute(nameof(Y)) != null) double.TryParse(xElement.Attribute(nameof(Y)).Value, NumberStyles.Any, CultureInfo.InvariantCulture, out y);
            X = x;
            Y = y;
            IsValid = true;
            if (double.IsInfinity(X) || double.IsNaN(X) || double.IsInfinity(Y) || double.IsNaN(Y)) { IsValid = false; }

        }
        #endregion

        #region Members

        /// <summary>
        /// X Value.
        /// </summary>
        public double X;

        /// <summary>
        /// Y Value.
        /// </summary>
        public double Y;

        /// <summary>
        /// Boolean indicating if the ordinate has valid numeric values or not.
        /// </summary>
        public bool IsValid;

        #endregion

        #region Methods
       
        /// <summary>
        /// Test if the ordinate is valid given monotonic criteria with the next/previous ordinate in a series.
        /// </summary>
        /// <param name="ordinateToCompare">Ordinate to compare to the target ordinate.</param>
        /// <param name="strictX">Are the x-values strictly monotonic?</param>
        /// <param name="strictY">Are the y-values strictly monotonic?</param>
        /// <param name="xOrder">The order of the x-values.</param>
        /// <param name="yOrder">The order of the y-values.</param>
        /// <param name="compareOrdinateIsNext">Boolean identifying if the ordinate to compare is the next or previous ordinate in a series.</param>
        /// <returns> A boolean indicating if the ordinate is valid or not given the criteria.</returns>
        public bool OrdinateValid(Ordinate ordinateToCompare, bool strictX, bool strictY, SortOrder xOrder, SortOrder yOrder, bool compareOrdinateIsNext)
        {
            // Check the ordinate itself
            if (IsValid == false) return false;
            // Check for monotonicity.
            if (compareOrdinateIsNext == true)
            {
                // Looking forward
                if (strictY && yOrder != SortOrder.None)
                {
                    if (ordinateToCompare.Y == Y)
                        return false;
                }

                if (yOrder == SortOrder.Descending)
                {
                    if (ordinateToCompare.Y > Y)
                        return false;
                }
                else if (yOrder == SortOrder.Ascending)
                {
                    if (ordinateToCompare.Y < Y)
                        return false;
                }

                if (strictX && xOrder != SortOrder.None)
                {
                    if (ordinateToCompare.X == X)
                        return false;
                }

                if (xOrder == SortOrder.Descending)
                {
                    if (ordinateToCompare.X > X)
                        return false;
                }
                else if (xOrder == SortOrder.Ascending)
                {
                    if (ordinateToCompare.X < X)
                        return false;
                }
            }
            else
            {
                // Looking back
                if (strictY && yOrder != SortOrder.None)
                {
                    if (Y == ordinateToCompare.Y)
                        return false;
                }

                if (yOrder == SortOrder.Descending)
                {
                    if (Y > ordinateToCompare.Y)
                        return false;
                }
                else if (yOrder == SortOrder.Ascending)
                {
                    if (Y < ordinateToCompare.Y)
                        return false;
                }

                if (strictX && xOrder != SortOrder.None)
                {
                    if (X == ordinateToCompare.X)
                        return false;
                }

                if (xOrder == SortOrder.Descending)
                {
                    if (X > ordinateToCompare.X)
                        return false;
                }
                else if (xOrder == SortOrder.Ascending)
                {
                    if (X < ordinateToCompare.X)
                        return false;
                }
            }
            // Passed the test
            return true;
        }

        /// <summary>
        /// Get any error messages with the ordinate given monotonic criteria with the next/previous ordinate in a series.
        /// </summary>
        /// <param name="ordinateToCompare">Ordinate to compare to the target ordinate.</param>
        /// <param name="strictX">Are the x-values strictly monotonic?</param>
        /// <param name="strictY">Are the y-values strictly monotonic?</param>
        /// <param name="xOrder">The order of the x-values.</param>
        /// <param name="yOrder">The order of the y-values.</param>
        /// <param name="compareOrdinateIsNext">Boolean identifying if the ordinate to compare is the next or previous ordinate in a series.</param>
        /// <returns> A list of error messages given the criteria.</returns>
        public List<string> OrdinateErrors(Ordinate ordinateToCompare, bool strictX, bool strictY, SortOrder xOrder, SortOrder yOrder, bool compareOrdinateIsNext)
        {
            var result = new List<string>();
            // Check the ordinate itself
            result.AddRange(OrdinateErrors());

            string orderY = (yOrder == SortOrder.Descending ? "decreasing" : "increasing");
            string orderX = (xOrder == SortOrder.Descending ? "decreasing" : "increasing");
            // Check for monotonicity.
            if (strictY && yOrder != SortOrder.None)
            {
                if (ordinateToCompare.Y == Y) { result.Add($"Y values must be strictly {orderY}."); }
            } 

            if (strictX && xOrder != SortOrder.None)
            {
                if (ordinateToCompare.X == X) { result.Add($"X values must be strictly {orderX}."); }
            } 

            if (compareOrdinateIsNext == true)
            {
                // Looking forward
                if (yOrder == SortOrder.Descending)
                {
                    if (ordinateToCompare.Y > Y) { result.Add("Y values must decrease."); }
                }
                else if (yOrder == SortOrder.Ascending)
                {
                    if (ordinateToCompare.Y < Y) { result.Add("Y values must increase."); }
                }

                if (xOrder == SortOrder.Descending)
                {
                    if (ordinateToCompare.X > X) { result.Add("X values must decrease."); }
                }
                else if (xOrder == SortOrder.Ascending)
                {
                    if (ordinateToCompare.X < X) { result.Add("X values must increase."); }
                }
            }
            else
            {
                // Looking back
                if (yOrder == SortOrder.Descending)
                {
                    if (Y > ordinateToCompare.Y) { result.Add("Y values must decrease."); }
                }
                else if (yOrder == SortOrder.Ascending)
                {
                    if (Y < ordinateToCompare.Y) { result.Add("Y values must increase."); }
                }

                if (xOrder == SortOrder.Descending)
                {
                    if (X > ordinateToCompare.X) { result.Add("X values must decrease."); }
                }
                else if (xOrder == SortOrder.Ascending)
                {
                    if (X < ordinateToCompare.X) { result.Add("X values must increase."); }
                }
            }
            // Done testing
            return result;
        }

        /// <summary>
        /// Get errors in the ordinate data.
        /// </summary>
        /// <returns>A list of strings.</returns>
        public List<string> OrdinateErrors()
        {
            var result = new List<string>();
            if (IsValid == false)
            {
                if (double.IsInfinity(X))
                    result.Add("X value can not be infinity.");
                if (double.IsNaN(X))
                    result.Add("X value must be a valid number.");
                if (double.IsInfinity(Y))
                    result.Add("Y value can not be infinity.");
                if (double.IsNaN(Y))
                    result.Add("Y value must be a valid number.");
            }
            return result;
        }

        /// <summary>
        /// Transform the x and y coordinates and return a transformed ordinate.
        /// </summary>
        /// <param name="xTransform">Transform method for the x value.</param>
        /// <param name="yTransform">Transform method for the y value.</param>
        public Ordinate Transform(Transform xTransform, Transform yTransform)
        {
            double transformedX = X, transformedY = Y;
            switch (xTransform)
            {
                case Data.Transform.Logarithmic:
                    {
                        transformedX = Tools.Log10(X);
                        break;
                    }
                case Data.Transform.NormalZ:
                    {
                        transformedX = Normal.StandardZ(X);
                        break;
                    }
            }
            switch (yTransform)
            {
                case Data.Transform.Logarithmic:
                    {
                        transformedY = Tools.Log10(Y);
                        break;
                    }
                case Data.Transform.NormalZ:
                    {
                        transformedY = Normal.StandardZ(Y);
                        break;
                    }
            }
            return new Ordinate(transformedX, transformedY);
        }

        /// <summary>
        /// Checks if two ordinates are equal.
        /// </summary>
        /// <param name="left"> First ordinate to compare.</param>
        /// <param name="right"> Second ordinate to compare.</param>
        /// <returns>True if two objects are numerically equal; otherwise, False.</returns>
        public static bool operator ==(Ordinate left, Ordinate right)
        {
            if (Math.Abs(left.X - right.X) > Tools.DoubleMachineEpsilon) { return false; }
            if (Math.Abs(left.Y - right.Y) > Tools.DoubleMachineEpsilon) { return false; }
            return true;
        }

        /// <summary>
        /// Checks if two ordinates are not equal.
        /// </summary>
        /// <param name="left"> First ordinate to compare.</param>
        /// <param name="right"> Second ordinate to compare.</param>
        /// <returns>True if two objects are not numerically equal; otherwise, False.</returns>
        public static bool operator !=(Ordinate left, Ordinate right)
        {
            return !(left == right);
        }

        /// <summary>
        /// Returns the ordinate as XElement.
        /// </summary>
        public XElement ToXElement()
        {
            var result = new XElement(nameof(Ordinate));
            result.SetAttributeValue(nameof(X), X.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(Y), Y.ToString("G17", CultureInfo.InvariantCulture));
            return result;
        }

        #endregion
    }
}