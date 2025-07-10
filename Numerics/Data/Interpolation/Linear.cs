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

using Numerics.Distributions;
using System;
using System.Collections.Generic;

namespace Numerics.Data
{
    /// <summary>
    /// A class for linear interpolation.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// Linear interpolation is a method of estimation and curve fitting using linear polynomials to construct new data points
    /// within the range of a discrete set of known data points.
    /// </para>
    /// <para>
    ///     Methods are designed to be compatible with the VBA Macro 'Interpolate Version 2.0.0' November 2017, USACE Risk Management Center
    /// </para>
    /// <b> References: </b>
    /// <see href="https://en.wikipedia.org/wiki/Linear_interpolation"/>
    /// </remarks>
    [Serializable]
    public class Linear : Interpolater
    {
        /// <summary>
        /// Construct new linear interpolation. 
        /// </summary>
        /// <param name="xValues">List of x-values.</param>
        /// <param name="yValues">List of y-values.</param>
        /// <param name="sortOrder">The sort order of the x-values, either ascending or descending. Default = Ascending. </param>
        public Linear(IList<double> xValues, IList<double> yValues, SortOrder sortOrder = SortOrder.Ascending) : base(xValues, yValues, sortOrder){}

        /// <summary>
        /// The transform for the x-values. Default = None.
        /// </summary>
        public Transform XTransform { get; set; } = Transform.None;

        /// <summary>
        /// The transform for the y-values. Default = None.
        /// </summary>
        public Transform YTransform { get; set; } = Transform.None;

        /// <inheritdoc/>
        public override double BaseInterpolate(double x, int index)
        {
            // See if x is out of range
            if ((SortOrder == SortOrder.Ascending && x <= XValues[0]) || 
                (SortOrder == SortOrder.Descending && x >= XValues[0])) return YValues[0];
            if ((SortOrder == SortOrder.Ascending && x >= XValues[Count - 1]) || 
                (SortOrder == SortOrder.Descending && x <= XValues[Count - 1])) return YValues[Count - 1];

            double y, x1 = 0, x2 = 0, y1 = 0, y2 = 0;
            int xlo = index, xhi = index + 1;

            // Get X transform
            if (XTransform == Transform.None)
            {
                x1 = XValues[xlo];
                x2 = XValues[xhi];
            }
            else if (XTransform == Transform.Logarithmic)
            {
                x = Tools.Log10(x);
                x1 = Tools.Log10(XValues[xlo]);
                x2 = Tools.Log10(XValues[xhi]);
            }
            else if (XTransform == Transform.NormalZ)
            {
                x = Normal.StandardZ(x);
                x1 = Normal.StandardZ(XValues[xlo]);
                x2 = Normal.StandardZ(XValues[xhi]);
            }

            // Get Y transform
            if (YTransform == Transform.None)
            {
                y1 = YValues[xlo];
                y2 = YValues[xhi];
            }
            else if (YTransform == Transform.Logarithmic)
            {
                y1 = Tools.Log10(YValues[xlo]);
                y2 = Tools.Log10(YValues[xhi]);
            }
            else if (YTransform == Transform.NormalZ)
            {
                y1 = Normal.StandardZ(YValues[xlo]);
                y2 = Normal.StandardZ(YValues[xhi]);
            }

            // Interpolate
            // Check for division by zero
            if ((x2 - x1) == 0)
            {
                y = y1;
            }
            else
            {
                y = y1 + (x - x1) / (x2 - x1) * (y2 - y1);
            }
            //
            if (YTransform == Transform.None)
            {
                return y;
            }
            else if (YTransform == Transform.Logarithmic)
            {
                return Math.Pow(10d, y);
            }
            else if (YTransform == Transform.NormalZ)
            {
                return Normal.StandardCDF(y);
            }

            // return NaN if we get to here
            return double.NaN;
        }

        /// <summary>
        /// Given a value x, returns an extrapolated value. 
        /// </summary>
        /// <param name="x">The value to extrapolate.</param>
        public double Extrapolate(double x)
        {
            double y, x1, x2, y1, y2;

            // See if x is out of range
            if ((SortOrder == SortOrder.Ascending && x < XValues[0]) || (SortOrder == SortOrder.Descending && x > XValues[0]))
            {
                x1 = XValues[0];
                x2 = XValues[1];
                y1 = YValues[0];
                y2 = YValues[1];
            }            
            else if ((SortOrder == SortOrder.Ascending && x > XValues[Count - 1]) || (SortOrder == SortOrder.Descending && x < XValues[Count - 1]))
            {
                x1 = XValues[Count - 1];
                x2 = XValues[Count - 2];
                y1 = YValues[Count - 1];
                y2 = YValues[Count - 2];
            }
            else
            {
                return Interpolate(x);
            }

            // Get X transform
            if (XTransform == Transform.Logarithmic)
            {
                x = Tools.Log10(x);
                x1 = Tools.Log10(x1);
                x2 = Tools.Log10(x2);
            }
            else if (XTransform == Transform.NormalZ)
            {
                x = Normal.StandardZ(x);
                x1 = Normal.StandardZ(x1);
                x2 = Normal.StandardZ(x2);
            }

            // Get Y transform
            if (YTransform == Transform.Logarithmic)
            {
                y1 = Tools.Log10(y1);
                y2 = Tools.Log10(y2);
            }
            else if (YTransform == Transform.NormalZ)
            {
                y1 = Normal.StandardZ(y1);
                y2 = Normal.StandardZ(y2);
            }

            // Extrapolate
            // Check for division by zero
            if ((y2 - y1) == 0)
            {
                y = y1;
            }
            else
            {
                y = y1 - (x1 - x) * (y2 - y1) / (x2 - x1);
            }
            //
            if (YTransform == Transform.None)
            {
                return y;
            }
            else if (YTransform == Transform.Logarithmic)
            {
                return Math.Pow(10d, y);
            }
            else if (YTransform == Transform.NormalZ)
            {
                return Normal.StandardCDF(y);
            }

            // return NaN if we get to here
            return double.NaN;
        }

    }

}
