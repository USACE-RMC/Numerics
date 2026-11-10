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
        /// Constructs a linear interpolator over the lists of another one, without validation.
        /// </summary>
        /// <param name="source">The interpolator to copy.</param>
        private Linear(Linear source) : base(source)
        {
            XTransform = source.XTransform;
            YTransform = source.YTransform;
        }

        /// <summary>
        /// Creates a copy that shares this interpolator's x and y lists.
        /// </summary>
        /// <returns>A new interpolator over the same lists, with the same sort order, transforms and search mode, and a fresh search state.</returns>
        /// <remarks>
        /// The lists are neither copied nor validated again, so the copy costs O(1). Interpolation only reads the lists,
        /// so each thread can use its own copy; a single instance must not be shared across threads, because every search
        /// updates its search state. The lists must not be modified while a copy is in use.
        /// </remarks>
        public Linear ShallowClone()
        {
            return new Linear(this);
        }

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
            // The NormalZ transform sends a value of 1 to +∞ (StandardZ(0) and Log10(0) stay finite).
            if (double.IsNaN(y)) y = OrderedPairedData.InfiniteEndpointLimit(x, x1, x2, y1, y2);
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

            // return NaN if this point is reached
            return double.NaN;
        }

        /// <summary>
        /// Given a list of x-values, returns an extrapolated value for each.
        /// </summary>
        /// <param name="x">The values to extrapolate.</param>
        /// <returns>One value per x-value, each equal to <see cref="Extrapolate(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="x"/> is null.</exception>
        /// <remarks>
        /// The values are evaluated in order through the scalar method, so the searches see the same sequence of
        /// points and every value has the same bits as the scalar call.
        /// </remarks>
        public double[] Extrapolate(IList<double> x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = Extrapolate(x[i]);
            return result;
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
            if ((x2 - x1) == 0)
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

            // return NaN if this point is reached
            return double.NaN;
        }

    }

}
