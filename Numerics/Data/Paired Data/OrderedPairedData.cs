﻿/**
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
* **/

using System;
using System.Collections;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Data;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Controls;
using System.Windows.Markup;
using System.Xml.Linq;
using Microsoft.VisualBasic;
using Numerics.Distributions;

namespace Numerics.Data
{

    /// <summary>
    /// The OrderedPairedData class is designed to store xy data that is ordered for both the x and y values.
    /// The most common use case for this class is to store curve data such as a cumulative distribution function.
    /// One of the more powerful functions associated with this class is the Transform() function which allows to transform
    /// the OrderedPairedData class with another to create a new class. Common use cases of the Transform() function include
    /// transforming a flood frequency function into a stage frequency function using a flow-stage rating curve.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     Authors:
    ///     Woodrow Fields, USACE Risk Management Center, woodrow.l.fields@usace.army.mil
    /// </para>
    /// </remarks> 
    [Serializable]
    public class OrderedPairedData : IList<Ordinate>, INotifyCollectionChanged
    {
        #region Search Properties

        /// <summary>
        /// Keeps track of the starting search location for x-values. 
        /// </summary>
        public int XSearchStart { get; set; }

        /// <summary>
        /// Keeps track of the starting search location for y-values. 
        /// </summary>
        public int YSearchStart { get; set; }

        /// <summary>
        /// Determines whether to use a smart searching algorithm or just sequential search.
        /// </summary>
        public bool UseSmartSearch { get; set; } = true;

        /// <summary>
        /// Keeps track of the difference is start locations. 
        /// </summary>
        private int XdeltaStart = 0;

        /// <summary>
        /// Keeps track of the difference is start locations. 
        /// </summary>
        private int YdeltaStart = 0;

        /// <summary>
        /// Determines which search method to use. If values are correlated, use the Hunt method. 
        /// </summary>
        private bool Xcorrelated = false;

        /// <summary>
        /// Determines which search method to use. If values are correlated, use the Hunt method. 
        /// </summary>
        private bool Ycorrelated = false;

        #endregion

        private bool _strictX;
        private bool _strictY;
        private SortOrder _orderX;
        private SortOrder _orderY;
        private List<Ordinate> _ordinates;

        public event NotifyCollectionChangedEventHandler CollectionChanged;

        /// <summary>
        /// Represents if the paired dataset has valid ordinates and order.
        /// </summary>
        /// <returns>boolean representing data validity.</returns>
        public bool IsValid { get; private set; }

        /// <summary>
        /// Number of ordinates in the dataset.
        /// </summary>
        /// <returns>Integer representing the number of ordinates.</returns>
        public int Count => _ordinates.Count;

        /// <summary>
        /// Determines whether the sort order is strict on the X variable. 
        /// </summary>
        public bool StrictX
        {
            get { return _strictX; }
            set
            {
                if (_strictX != value)
                {
                    _strictX = value;
                    Validate();
                }
            }
        }

        /// <summary>
        /// Determines whether the sort order is strict on the Y variable. 
        /// </summary>
        public bool StrictY
        {
            get { return _strictY; }
            set
            {
                if (_strictY != value)
                {
                    _strictY = value;
                    Validate();
                }
            }
        }

        /// <summary>
        /// Gets or sets the sort order of the X variable. 
        /// </summary>
        public SortOrder OrderX
        {
            get { return _orderX; }
            set
            {
                if (_orderX != value)
                {
                    _orderX = value;
                    Validate();
                }
            }
        }

        /// <summary>
        /// Gets or sets the sort order of the Y variable. 
        /// </summary>
        public SortOrder OrderY
        {
            get { return _orderY; }
            set
            {
                if (_orderY != value)
                {
                    _orderY = value;
                    Validate();
                }
            }
        }

        /// <summary>
        /// ReadOnly is an implementation of ICollection and not implemented for this class.
        /// </summary>
        /// <returns>Boolean.</returns>
        bool ICollection<Ordinate>.IsReadOnly => false;

        /// <summary>
        /// Get or set the ordinate at a specified index.
        /// </summary>
        /// <param name="index">Index of ordinate in the collection to manipulate.</param>
        /// <returns>Ordinate at the specified index.</returns>
        public Ordinate this[int index]
        {
            get { return _ordinates[index]; }
            set
            {
                //_ordinates[index] = value;
                //Validate();

                if (_ordinates[index] != value)
                {
                    var oldvalue = _ordinates[index];
                    _ordinates[index] = value;
                    if (IsValid == true)
                    {
                        if (OrdinateValid(index) == false) { IsValid = false; }
                    }
                    else
                    {
                        if (OrdinateValid(index) == true) { Validate(); }
                    }
                    ////
                    //if (SupressCollectionChanged == false)
                    //{
                    CollectionChanged?.Invoke(this, new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Replace, value, oldvalue));
                    //}
                }

            }
        }

        /// <summary>
        /// Create empty instance of the ordered paired data class.
        /// </summary>
        /// <param name="strictOnX">Boolean indicating if x values are strictly increasing/decreasing. True means x values cannot be equal.</param>
        /// <param name="xOrder">Order of the x values.</param>
        /// <param name="strictOnY">Boolean indicating if x values are strictly increasing/decreasing. True means x values cannot be equal.</param>
        /// <param name="yOrder">Order of the y values.</param>
        public OrderedPairedData(bool strictOnX, SortOrder xOrder, bool strictOnY, SortOrder yOrder)
        {
            _ordinates = new List<Ordinate>();
            StrictX = strictOnX;
            StrictY = strictOnY;
            OrderX = xOrder;
            OrderY = yOrder;
        }

        /// <summary>
        /// Create empty instance of the ordered paired data class with a preset capacity.
        /// </summary>
        /// <param name="capacity">capacity of the collection.</param>
        /// <param name="strictOnX">Boolean indicating if x values are strictly increasing/decreasing. True means x values cannot be equal.</param>
        /// <param name="xOrder">Order of the x values.</param>
        /// <param name="strictOnY">Boolean indicating if x values are strictly increasing/decreasing. True means x values cannot be equal.</param>
        /// <param name="yOrder">Order of the y values.</param>
        public OrderedPairedData(int capacity, bool strictOnX, SortOrder xOrder, bool strictOnY, SortOrder yOrder)
        {
            _ordinates = new List<Ordinate>(capacity);
            StrictX = strictOnX;
            StrictY = strictOnY;
            OrderX = xOrder;
            OrderY = yOrder;
        }

        /// <summary>
        /// Create an instance of the ordered paired data class with defined ordinate data.
        /// </summary>
        /// <param name="xData">Ordinate x values.</param>
        /// <param name="yData">Ordinate y values.</param>
        /// <param name="strictOnX">Boolean indicating if x values are strictly increasing/decreasing. True means x values cannot be equal.</param>
        /// <param name="xOrder">Order of the x values.</param>
        /// <param name="strictOnY">Boolean indicating if x values are strictly increasing/decreasing. True means x values cannot be equal.</param>
        /// <param name="yOrder">Order of the y values.</param>
        public OrderedPairedData(IList<double> xData, IList<double> yData, bool strictOnX, SortOrder xOrder, bool strictOnY, SortOrder yOrder)
        {
            StrictX = strictOnX;
            StrictY = strictOnY;
            OrderX = xOrder;
            OrderY = yOrder;
            _ordinates = new List<Ordinate>(xData.Count);
            for (int i = 0; i < xData.Count; i++)
                _ordinates.Add(new Ordinate(xData[i], yData[i]));
            Validate();
        }

        /// <summary>
        /// Create an instance of the ordered paired data class with defined ordinate data.
        /// </summary>
        /// <param name="data">Ordinate values.</param>
        /// <param name="strictOnX">Boolean indicating if x values are strictly increasing/decreasing. True means x values cannot be equal.</param>
        /// <param name="xOrder">Order of the x values.</param>
        /// <param name="strictOnY">Boolean indicating if x values are strictly increasing/decreasing. True means x values cannot be equal.</param>
        /// <param name="yOrder">Order of the y values.</param>
        public OrderedPairedData(IList<Ordinate> data, bool strictOnX, SortOrder xOrder, bool strictOnY, SortOrder yOrder)
        {
            StrictX = strictOnX;
            StrictY = strictOnY;
            OrderX = xOrder;
            OrderY = yOrder;
            _ordinates = new List<Ordinate>(data.Count);
            for (int i = 0; i < data.Count; i++)
                _ordinates.Add(new Ordinate(data[i].X, data[i].Y));
            Validate();
        }

        /// <summary>
        /// Create a new instance of the uncertain ordered paired data class from an XElement XML object.
        /// </summary>
        /// <param name="el"></param>
        public OrderedPairedData(XElement el)
        {
            // Get Strictness
            bool strict = false;
            if (el.Attribute(nameof(StrictX)) != null) { bool.TryParse(el.Attribute(nameof(StrictX)).Value, out strict); }
            StrictX = strict;

            strict = false;
            if (el.Attribute(nameof(StrictY)) != null) { bool.TryParse(el.Attribute(nameof(StrictY)).Value, out strict); }
            StrictY = strict;

            // Get Order
            SortOrder order = SortOrder.None;
            if (el.Attribute(nameof(OrderX)) != null) { Enum.TryParse(el.Attribute(nameof(OrderX)).Value, out order); }
            OrderX = order;

            order = SortOrder.None;
            if (el.Attribute(nameof(OrderY)) != null) { Enum.TryParse(el.Attribute(nameof(OrderY)).Value, out order); }
            OrderY = order;

            // Ordinates
            var curveEl = el.Element("Ordinates");
            _ordinates = new List<Ordinate>();
            if (curveEl != null)
            {
                foreach (XElement o in curveEl.Elements(nameof(Ordinate)))
                {
                    _ordinates.Add(new Ordinate(o));
                }
            }

            Validate();
        }



        /// <summary>
        /// Determines the index of a specific ordinate in the collection.
        /// </summary>
        /// <param name="item">The ordinate to locate in the collection.</param>
        /// <returns>The index of ordinate item if found in the list; otherwise, -1.</returns>
        public int IndexOf(Ordinate item)
        {
            return _ordinates.IndexOf(item);
        }

        /// <summary>
        /// Determines the first index of an ordinate where the x and y values are equal in the collection.
        /// </summary>
        /// <param name="xValue">The X value to locate in the collection.</param>
        /// <param name="yValue">The Y value to locate in the collection.</param>
        /// <returns>The index of the first ordinate item with equivalent x and y values if found in the list; otherwise, -1.</returns>
        public int IndexOf(double xValue, double yValue)
        {
            // This function could be really good for getting the index based on the values of the ordinates. It could burn if there are multiple items with the same values and someone is trying to remove a specific ordinate.
            for (int i = 0; i < _ordinates.Count; i++)
            {
                if (Math.Abs(_ordinates[i].X - xValue) <= 0.000000001d && Math.Abs(_ordinates[i].Y - yValue) <= 0.000000001d) return i;
            }
            return -1;
        }

        /// <summary>
        /// Method to determine validity of all ordinates in the collection.
        /// </summary>
        private void Validate()
        {
            if (_ordinates == null) return;
            IsValid = true;
            for (int i = 0; i < _ordinates.Count; i++)
            {
                if (OrdinateValid(i, false) == false)
                {
                    IsValid = false;
                    return;
                }
            }
        }

        /// <summary>
        /// Determines the validity of a specific ordinate at a specified index.
        /// </summary>
        /// <param name="index">The index of ordinate item.</param>
        /// <param name="lookBackward">Optional parameter to also look backward when determining validity of the ordinate. This parameter is included as an optimization for situations where looking at the previous ordinate is not required.</param>
        /// <returns>A boolean that indicates if the ordinate at the specified index is valid or not within the dataset.</returns>
        private bool OrdinateValid(int index, bool lookBackward = true)
        {
            if (index < 0 || index > _ordinates.Count - 1)
                return true;
            // Look backward
            if (lookBackward == true & index > 0)
            {
                if (_ordinates[index].OrdinateValid(_ordinates[index - 1], StrictX, StrictY, OrderX, OrderY, false) == false)
                    return false;
            }
            // Look forward
            if (index < _ordinates.Count - 1)
            {
                if (_ordinates[index].OrdinateValid(_ordinates[index + 1], StrictX, StrictY, OrderX, OrderY, true) == false)
                    return false;
            }
            // Passed the test
            return true;
        }

        /// <summary>
        /// Get a list of errors in the UncertainOrderedPairedData object.
        /// </summary>
        /// <returns>A list of strings.</returns>
        public List<string> GetErrors()
        {
            var result = new List<string>();
            if (IsValid) return result;
            for (int i = 0; i < _ordinates.Count - 1; i++)
                // Look forward
                result.AddRange(_ordinates[i].OrdinateErrors(_ordinates[i + 1], StrictX, StrictY, OrderX, OrderY, true));
            result.AddRange(_ordinates.Last().OrdinateErrors());
            return result;
        }

        /// <summary>
        /// Removes the first occurrence of a target ordinate from the collection.
        /// </summary>
        /// <param name="item">The ordinate to remove from the collection.</param>
        /// <returns>True if item was successfully removed from the collection; otherwise, False.</returns>
        public bool Remove(Ordinate item)
        {
            int itemIndex = IndexOf(item);
            if (itemIndex == -1) return false;
            _ordinates.RemoveAt(itemIndex);
            Validate();
            CollectionChanged?.Invoke(this, new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Remove, item, itemIndex));
            return true;
        }

        /// <summary>
        /// Removes the ordinate from the collection at the specified index.
        /// </summary>
        /// <param name="index">The zero-based index of the item to remove.</param>
        public void RemoveAt(int index)
        {
            if (index < 0 || index >= _ordinates.Count) { return; }
            var item = _ordinates[index];
            _ordinates.RemoveAt(index);
            Validate();
            CollectionChanged?.Invoke(this, new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Remove, item, index));
        }

        /// <summary>
        /// Removes a range of ordinates from the collection.
        /// </summary>
        /// <param name="index">The zero-based starting index of the range of elements to remove.</param>
        /// <param name="count">The number of elements to remove.</param>
        public void RemoveRange(int index, int count)
        {
            if (index < 0 || (index + count) >= _ordinates.Count) { return; }
            List<Ordinate> items = new List<Ordinate>();
            for (int i = index; i < count; i++) { items.Add(_ordinates[i]); }
            _ordinates.RemoveRange(index, count);
            Validate();
            CollectionChanged?.Invoke(this, new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Remove, items, index));
        }

        /// <summary>
        /// Adds an ordinate to the collection.
        /// </summary>
        /// <param name="item">The ordinate to add to the collection.</param>
        public void Add(Ordinate item)
        {
            _ordinates.Add(item);
            IsValid = OrdinateValid(_ordinates.Count - 1);
            CollectionChanged?.Invoke(this, new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Add, item, _ordinates.Count - 1));
        }

        /// <summary>
        /// Inserts an item to the collection of ordinates.
        /// </summary>
        /// <param name="index">The zero-based index at which item should be inserted.</param>
        /// <param name="item">The ordinate to insert into the collection.</param>
        public void Insert(int index, Ordinate item)
        {
            _ordinates.Insert(index, item);
            // only need to set valid state if it is true..if it is already false then inserting can't make it true.
            if (IsValid) IsValid = OrdinateValid(index);
            CollectionChanged?.Invoke(this, new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Add, item, index));
        }

        /// <summary>
        /// Removes all ordinates form the collection and sets count to zero.
        /// </summary>
        public void Clear()
        {
            _ordinates.Clear();
            IsValid = true;
            CollectionChanged?.Invoke(this, new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Reset));
        }

        /// <summary>
        /// Determines whether the collection contains a specific ordinate.
        /// </summary>
        /// <param name="item">The ordinate to locate in the collection.</param>
        /// <returns>True if ordinate is found in the collection; otherwise, false.</returns>
        public bool Contains(Ordinate item)
        {
            return IndexOf(item) != -1;
        }

        /// <summary>
        /// Calculates the area between the Y values and the zero axis using the trapezoidal approximation.
        /// </summary>
        /// <returns>The area under the curve.</returns>
        public double TrapezoidalAreaUnderY()
        {
            double sum = 0d;
            if (OrderX == SortOrder.Ascending)
            {
                for (int i = 1; i < _ordinates.Count; i++)
                    sum += 0.5d * (_ordinates[i].X - _ordinates[i - 1].X) * (_ordinates[i - 1].Y + _ordinates[i].Y);
            }
            else if (OrderX == SortOrder.Descending)
            {
                for (int i = 1; i < _ordinates.Count; i++)
                    sum += 0.5d * (_ordinates[i - 1].X - _ordinates[i].X) * (_ordinates[i - 1].Y + _ordinates[i].Y);
            }
            else
            {
                throw new InvalidConstraintException("Unable to calculate area under y with no guarantee of sort order on x.");
            }
            return sum;
        }

        /// <summary>
        /// Calculates the area between the X values and the zero axis using the trapezoidal approximation.
        /// </summary>
        /// <returns>The area under the curve.</returns>
        public double TrapezoidalAreaUnderX()
        {
            double sum = 0d;
            if (OrderY == SortOrder.Ascending)
            {
                for (int i = 1; i < _ordinates.Count; i++)
                    sum += 0.5d * (_ordinates[i].Y - _ordinates[i - 1].Y) * (_ordinates[i - 1].X + _ordinates[i].X);
            }
            else if (OrderY == SortOrder.Descending)
            {
                for (int i = 1; i < _ordinates.Count; i++)
                    sum += 0.5d * (_ordinates[i - 1].Y - _ordinates[i].Y) * (_ordinates[i - 1].X + _ordinates[i].X);
            }
            else
            {
                throw new InvalidConstraintException("Unable to calculate area under x with no guarantee of sort order on y.");
            }
            return sum;
        }

        #region Interpolation


        private double RawInterpolate(double value, int start, bool givenX = true, Transform xTransform = Data.Transform.None, Transform yTransform = Data.Transform.None)
        {

            double x = givenX ? value : 0, y = givenX ? 0 : value, x1 = 0, x2 = 0, y1 = 0, y2 = 0;
            int lo = start, hi = start + 1;

            // Get X transform
            if (xTransform == Data.Transform.None)
            {
                x1 = _ordinates[lo].X;
                x2 = _ordinates[hi].X;
            }
            else if (xTransform == Data.Transform.Logarithmic)
            {
                x = Tools.Log10(x);
                x1 = Tools.Log10(_ordinates[lo].X);
                x2 = Tools.Log10(_ordinates[hi].X);
            }
            else if (xTransform == Data.Transform.NormalZ)
            {
                x = Normal.StandardZ(x);
                x1 = Normal.StandardZ(_ordinates[lo].X);
                x2 = Normal.StandardZ(_ordinates[hi].X);
            }

            // Get Y transform
            if (yTransform == Data.Transform.None)
            {
                y1 = _ordinates[lo].Y;
                y2 = _ordinates[hi].Y;
            }
            else if (yTransform == Data.Transform.Logarithmic)
            {
                y = Tools.Log10(y);
                y1 = Tools.Log10(_ordinates[lo].Y);
                y2 = Tools.Log10(_ordinates[hi].Y);
            }
            else if (yTransform == Data.Transform.NormalZ)
            {
                y = Normal.StandardZ(y);
                y1 = Normal.StandardZ(_ordinates[lo].Y);
                y2 = Normal.StandardZ(_ordinates[hi].Y);
            }

            // interpolate
            if (givenX)
            {
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
                if (yTransform == Data.Transform.None)
                {
                    return y;
                }
                else if (yTransform == Data.Transform.Logarithmic)
                {
                    return Math.Pow(10d, y);
                }
                else if (yTransform == Data.Transform.NormalZ)
                {
                    return Normal.StandardCDF(y);
                }
            }
            else
            {
                // Check for division by zero
                if ((y2 - y1) == 0)
                {
                    x = x1;
                }
                else
                {
                    x = x1 + (y - y1) / (y2 - y1) * (x2 - x1);
                }
                //
                if (xTransform == Data.Transform.None)
                {
                    return x;
                }
                else if (xTransform == Data.Transform.Logarithmic)
                {
                    return Math.Pow(10d, x);
                }
                else if (xTransform == Data.Transform.NormalZ)
                {
                    return Normal.StandardCDF(x);
                }
            }

            // return NaN if we get to here
            return double.NaN;

        }

        public double Interpolate(double value, bool givenX = true, Transform xTransform = Data.Transform.None, Transform yTransform = Data.Transform.None)
        {
            if (Count == 0) return double.NaN;

            // First see if value is out of range
            if (givenX)
            {
                if (OrderX == SortOrder.None)
                    throw new Exception("Interpolation requires the x-values to be ascending or descending.");

                if (Count == 1) return _ordinates[0].Y;
                if ((OrderX == SortOrder.Ascending && value <= _ordinates[0].X) || (OrderX == SortOrder.Descending && value >= _ordinates[0].X)) return _ordinates[0].Y;
                if ((OrderX == SortOrder.Ascending && value >= _ordinates[Count - 1].X) || (OrderX == SortOrder.Descending && value <= _ordinates[Count - 1].X)) return _ordinates[Count - 1].Y;
            }
            else
            {
                if (OrderY == SortOrder.None)
                    throw new Exception("Interpolation requires the y-values to be ascending or descending.");

                if (Count == 1) return _ordinates[0].X;
                if ((OrderY == SortOrder.Ascending && value <= _ordinates[0].Y) || (OrderY == SortOrder.Descending && value >= _ordinates[0].Y)) return _ordinates[0].X;
                if ((OrderY == SortOrder.Ascending && value >= _ordinates[Count - 1].Y) || (OrderY == SortOrder.Descending && value <= _ordinates[Count - 1].Y)) return _ordinates[Count - 1].X;
            }
            // Interpolate
            int start = givenX ? SearchX(value) : SearchY(value);
            return RawInterpolate(value, start, givenX, xTransform, yTransform);
        }

        public double[] Interpolate(IList<double> values, bool givenX = true, Transform xTransform = Data.Transform.None, Transform yTransform = Data.Transform.None)
        {
            var result = new double[values.Count];
            for (int i = 0; i < values.Count; i++)
                result[i] = Interpolate(values[i], givenX, xTransform, yTransform);
            return result;
        }


        /// <summary>
        /// Samples a Y value for a given X from the X coordinates of the curve, solves the function F(X) = Y
        /// </summary>
        /// <param name="xValue">A value that represents the X, if the value is below the lowest x value, it returns the first y value, if the value is above the highest x value it returns the last y value</param>
        /// <param name="xTransform">Interpolation transform for x values.</param>
        /// <param name="yTransform">Interpolation transform for y values.</param>
        /// <returns>An interpolated y value for a given x</returns>
        public double GetYfromX(double xValue, Transform xTransform = Data.Transform.None, Transform yTransform = Data.Transform.None)
        {
            int index = BinarySearchX(xValue);
            if (index < -1)
            {
                index = -1 * index - 1;
                if (index == _ordinates.Count)
                {
                    return _ordinates.Last().Y;
                }
                else
                {
                    return InterpolateYFromX(xValue, _ordinates[index - 1], _ordinates[index], xTransform, yTransform);
                }
            }
            else if (index == -1)
            {
                return _ordinates[0].Y;
            }
            else
            {
                return _ordinates[index].Y;
            }
        }

        /// <summary>
        /// Interpolate Y values for a collection of x values.
        /// </summary>
        /// <param name="xValues">IList of x values to interpolate y values from.</param>
        /// <param name="valuesOrdered">Boolean indicating if the sample x values are in the same order as the ordinate x values.</param>
        /// <param name="xTransform">Interpolation transform for x values.</param>
        /// <param name="yTransform">Interpolation transform for y values.</param>
        public double[] GetYValues(IList<double> xValues, bool valuesOrdered, Transform xTransform = Data.Transform.None, Transform yTransform = Data.Transform.None)
        {
            if (xValues == null || xValues.Count == 0) return new double[0];
            double[] result = new double[xValues.Count];

            if (valuesOrdered == false)
            {
                //Sort values and keep track of original indices
                var sorted = xValues.Select((x, i) => new KeyValuePair<double, int>(x, i)).OrderBy(x => x.Key).ToList();
                if (OrderX == SortOrder.Descending) sorted.Reverse();
                //Get sorted values and original index locations
                double[] B = sorted.Select(x => x.Key).ToArray();
                int[] idx = sorted.Select(x => x.Value).ToArray();
                //Sample y values with sorted x and return in the original order.
                double[] sortedResults = GetYValues(B, true, xTransform, yTransform);
                for (int i = 0; i < sortedResults.Length; i++)
                    result[idx[i]] = sortedResults[i];
                return result;
            }

            // Get initial search index from first value
            int searchIndex = BinarySearchX(xValues[0]);
            if (searchIndex <= -1) searchIndex = -1 * searchIndex - 1;
            result[0] = (searchIndex == 0) ? _ordinates[0].Y : InterpolateYFromX(xValues[0], _ordinates[searchIndex - 1], _ordinates[searchIndex], xTransform, yTransform);
            //
            double xValue;
            if (OrderX == SortOrder.Ascending)
            {
                for (int i = 1; i < xValues.Count; i++)
                {
                    xValue = xValues[i];
                    for (int j = searchIndex; j < _ordinates.Count; j++)
                    {
                        if (_ordinates[j].X >= xValue) { searchIndex = j; break; }
                    }
                    result[i] = (searchIndex == 0) ? _ordinates[0].Y : InterpolateYFromX(xValue, _ordinates[searchIndex - 1], _ordinates[searchIndex], xTransform, yTransform);
                }
            }
            else
            {
                for (int i = 1; i < xValues.Count; i++)
                {
                    xValue = xValues[i];
                    for (int j = searchIndex; j < _ordinates.Count; j++)
                    {
                        if (_ordinates[j].X <= xValue) { searchIndex = j; break; }
                    }
                    result[i] = (searchIndex == 0) ? _ordinates[0].Y : InterpolateYFromX(xValue, _ordinates[searchIndex - 1], _ordinates[searchIndex], xTransform, yTransform);
                }
            }

            return result;
        }

        /// <summary>
        /// Samples a X value based on a given Y from the Y coordinates of the curve, solves the inverse function of F(X)=Y
        /// </summary>
        /// <param name="yValue">A value that represents the Y, if the value is below the lowest y value, it returns the lowest x value, if the value is above the highest y value it returns the highest x value</param>
        /// <param name="xTransform">Interpolation transform for x values.</param>
        /// <param name="yTransform">Interpolation transform for y values.</param>
        /// <returns>an x value for a given y</returns>
        /// <remarks></remarks>
        public double GetXfromY(double yValue, Transform xTransform = Data.Transform.None, Transform yTransform = Data.Transform.None)
        {
            int index = BinarySearchY(yValue);
            if (index < -1)
            {
                index = -1 * index - 1;
                if (index == _ordinates.Count)
                {
                    return _ordinates.Last().X;
                }
                else
                {
                    return InterpolateXFromY(yValue, _ordinates[index - 1], _ordinates[index], xTransform, yTransform);
                }
            }
            else if (index == -1)
            {
                return _ordinates[0].X;
            }
            else
            {
                return _ordinates[index].X;
            }
        }





        /// <summary>
        /// Calculates the interpolated Y value for a given X between two ordinates.
        /// </summary>
        /// <param name="x">Given X value to interpolate between two ordinates.</param>
        /// <param name="p1">First ordinate point.</param>
        /// <param name="p2">Second ordinate point.</param>
        /// <param name="xTransform">Interpolation transform for x values.</param>
        /// <param name="yTransform">Interpolation transform for y values.</param>
        /// <returns>Y value for the given X between two points.</returns>
        public static double InterpolateYFromX(double x, Ordinate p1, Ordinate p2, Transform xTransform, Transform yTransform)
        {
            double xValue = default, x1 = default, x2 = default;
            double y1, y2;
            // 
            switch (xTransform)
            {
                case Data.Transform.None:
                    {
                        xValue = x;
                        x1 = p1.X;
                        x2 = p2.X;
                        break;
                    }

                case Data.Transform.Logarithmic:
                    {
                        xValue = Tools.Log10(x);
                        x1 = Tools.Log10(p1.X);
                        x2 = Tools.Log10(p2.X);
                        break;
                    }

                case Data.Transform.NormalZ:
                    {
                        xValue = Normal.StandardZ(x);
                        x1 = Normal.StandardZ(p1.X);
                        x2 = Normal.StandardZ(p2.X);
                        break;
                    }
            }
            // 
            switch (yTransform)
            {
                case Data.Transform.None:
                    {
                        y1 = p1.Y;
                        y2 = p2.Y;
                        return y1 + (xValue - x1) / (x2 - x1) * (y2 - y1);
                    }

                case Data.Transform.Logarithmic:
                    {
                        y1 = Tools.Log10(p1.Y);
                        y2 = Tools.Log10(p2.Y);
                        return Math.Pow(10d, y1 + (xValue - x1) / (x2 - x1) * (y2 - y1));
                    }

                case Data.Transform.NormalZ:
                    {
                        y1 = Normal.StandardZ(p1.Y);
                        y2 = Normal.StandardZ(p2.Y);
                        return Normal.StandardCDF(y1 + (xValue - x1) / (x2 - x1) * (y2 - y1));
                    }
            }

            return default;
        }

        /// <summary>
        /// Calculates the interpolated X value for a given Y between two ordinates.
        /// </summary>
        /// <param name="y">Given Y value to interpolate between two ordinates.</param>
        /// <param name="p1">First ordinate point.</param>
        /// <param name="p2">Second ordinate point.</param>
        /// <param name="xTransform">Interpolation transform for x values.</param>
        /// <param name="yTransform">Interpolation transform for y values.</param>
        /// <returns>X value for the given Y between two points.</returns>
        public static double InterpolateXFromY(double y, Ordinate p1, Ordinate p2, Transform xTransform, Transform yTransform)
        {
            double x1 = default, x2 = default;
            double yValue, y1, y2;
            // 
            switch (xTransform)
            {
                case Data.Transform.None:
                    {
                        x1 = p1.X;
                        x2 = p2.X;
                        break;
                    }

                case Data.Transform.Logarithmic:
                    {
                        x1 = Tools.Log10(p1.X);
                        x2 = Tools.Log10(p2.X);
                        break;
                    }

                case Data.Transform.NormalZ:
                    {
                        x1 = Normal.StandardZ(p1.X);
                        x2 = Normal.StandardZ(p2.X);
                        break;
                    }
            }
            // 
            switch (yTransform)
            {
                case Data.Transform.None:
                    {
                        yValue = y;
                        y1 = p1.Y;
                        y2 = p2.Y;
                        return x1 + (yValue - y1) / (y2 - y1) * (x2 - x1);
                    }

                case Data.Transform.Logarithmic:
                    {
                        yValue = Tools.Log10(y);
                        y1 = Tools.Log10(p1.Y);
                        y2 = Tools.Log10(p2.Y);
                        return Math.Pow(10d, x1 + (yValue - y1) / (y2 - y1) * (x2 - x1));
                    }

                case Data.Transform.NormalZ:
                    {
                        yValue = Normal.StandardZ(y);
                        y1 = Normal.StandardZ(p1.Y);
                        y2 = Normal.StandardZ(p2.Y);
                        return Normal.StandardCDF(x1 + (yValue - y1) / (y2 - y1) * (x2 - x1));
                    }
            }

            return default;
            // 
            // Return x1 + (yValue - y1) / (y2 - y1) * (x2 - x1)
        }

        /// <summary>
        /// Transforms target with another OrderedPairedData collection into a new OrderedPairedData collection.
        /// </summary>
        /// <param name="transformFunction">The OrderedPairedData collection to be used for composition with target.</param>
        /// <param name="sourceCommonVariableX">Defines which ordinate value (x or y) the target function shares with the transform.</param>
        /// <param name="transformCommonVariableX">Defines which ordinate value (x or y) the transform function shares with the target.</param>
        /// <param name="xTransform">Interpolation transform for x values.</param>
        /// <param name="yTransform">Interpolation transform for y values.</param>
        public OrderedPairedData Transform(OrderedPairedData transformFunction, bool sourceCommonVariableX = true, bool transformCommonVariableX = true, Transform xTransform = Data.Transform.None, Transform yTransform = Data.Transform.None)
        {
            // Determine the starting index and the step direction which changes depending on ascending or descending order.
            int sourceStep = 1;
            int sourceIndex = 0;
            if (sourceCommonVariableX == true & OrderX == SortOrder.None)
                throw new ConstraintException("x values of source must be in a sorted order to transform.");
            if (transformCommonVariableX == true & transformFunction.OrderX == SortOrder.None)
                throw new ConstraintException("x values of transform must be in a sorted order to transform.");
            if (sourceCommonVariableX == false & OrderY == SortOrder.None)
                throw new ConstraintException("y values of source must be in a sorted order to transform.");
            if (transformCommonVariableX == false & transformFunction.OrderY == SortOrder.None)
                throw new ConstraintException("y values of transform must be in a sorted order to transform.");
            if (sourceCommonVariableX == true & OrderX == SortOrder.Descending)
            {
                sourceStep = -1;
                sourceIndex = _ordinates.Count - 1;
            }

            if (sourceCommonVariableX == false & OrderY == SortOrder.Descending)
            {
                sourceStep = -1;
                sourceIndex = _ordinates.Count - 1;
            }
            // 
            var tOrdinates = transformFunction._ordinates;
            int transformStep = 1;
            int transformIndex = 0;
            if (transformCommonVariableX == true & transformFunction.OrderX == SortOrder.Descending)
            {
                transformStep = -1;
                transformIndex = tOrdinates.Count - 1;
            }

            if (transformCommonVariableX == false & transformFunction.OrderY == SortOrder.Descending)
            {
                transformStep = -1;
                transformIndex = tOrdinates.Count - 1;
            }
            // 
            // Find the first index to start at, need to ignore non-overlapping sections.

            if (sourceCommonVariableX & transformCommonVariableX)
            {
                if (_ordinates[sourceIndex].X > tOrdinates[transformIndex].X)
                {
                    transformIndex += transformStep;
                    while (_ordinates[sourceIndex].X >= tOrdinates[transformIndex].X)
                        transformIndex += transformStep;
                }
                else
                {
                    sourceIndex += sourceStep;
                    while (_ordinates[sourceIndex].X <= tOrdinates[transformIndex].X)
                        sourceIndex += sourceStep;
                }
            }
            else if (sourceCommonVariableX & transformCommonVariableX == false)
            {
                if (_ordinates[sourceIndex].X > tOrdinates[transformIndex].Y)
                {
                    transformIndex += transformStep;
                    while (_ordinates[sourceIndex].X >= tOrdinates[transformIndex].Y)
                        transformIndex += transformStep;
                }
                else
                {
                    sourceIndex += sourceStep;
                    while (_ordinates[sourceIndex].X <= tOrdinates[transformIndex].Y)
                        sourceIndex += sourceStep;
                }
            }
            else if (sourceCommonVariableX == false & transformCommonVariableX)
            {
                if (_ordinates[sourceIndex].Y > tOrdinates[transformIndex].X)
                {
                    transformIndex += transformStep;
                    while (_ordinates[sourceIndex].Y >= tOrdinates[transformIndex].X)
                        transformIndex += transformStep;
                }
                else
                {
                    sourceIndex += sourceStep;
                    while (_ordinates[sourceIndex].Y <= tOrdinates[transformIndex].X)
                        sourceIndex += sourceStep;
                }
            }
            else if (sourceCommonVariableX == false & transformCommonVariableX == false)
            {
                if (_ordinates[sourceIndex].Y > tOrdinates[transformIndex].Y)
                {
                    transformIndex += transformStep;
                    while (_ordinates[sourceIndex].Y >= tOrdinates[transformIndex].Y)
                        transformIndex += transformStep;
                }
                else
                {
                    sourceIndex += sourceStep;
                    while (_ordinates[sourceIndex].Y <= tOrdinates[transformIndex].Y)
                        sourceIndex += sourceStep;
                }
            }
            // compose the two curves up to the overlapping areas.
            var xValues = new List<double>(); // always from source
            var yValues = new List<double>(); // always from target
            if (sourceCommonVariableX & transformCommonVariableX)
            {
                while (!(sourceIndex < 0 | sourceIndex == _ordinates.Count | transformIndex < 0 | transformIndex == tOrdinates.Count))
                {
                    if (_ordinates[sourceIndex].X == tOrdinates[transformIndex].X)
                    {
                        xValues.Add(_ordinates[sourceIndex].Y);
                        yValues.Add(tOrdinates[transformIndex].Y);
                        transformIndex += transformStep;
                        sourceIndex += sourceStep;
                    }
                    else if (_ordinates[sourceIndex].X > tOrdinates[transformIndex].X)
                    {
                        xValues.Add(InterpolateYFromX(tOrdinates[transformIndex].X, _ordinates[sourceIndex - 1], _ordinates[sourceIndex], xTransform, yTransform)); // GetYfromX(tOrdinates(transformIndex).X))
                        yValues.Add(tOrdinates[transformIndex].Y);
                        transformIndex += transformStep;
                    }
                    else
                    {
                        xValues.Add(_ordinates[sourceIndex].Y);
                        yValues.Add(InterpolateYFromX(_ordinates[sourceIndex].X, tOrdinates[transformIndex - 1], tOrdinates[transformIndex], xTransform, yTransform)); // transform.GetYfromX(_ordinates(sourceIndex).X))
                        sourceIndex += sourceStep;
                    }
                }
            }
            else if (sourceCommonVariableX & transformCommonVariableX == false)
            {
                while (!(sourceIndex < 0 | sourceIndex == _ordinates.Count | transformIndex < 0 | transformIndex == tOrdinates.Count))
                {
                    if (_ordinates[sourceIndex].X == tOrdinates[transformIndex].Y)
                    {
                        xValues.Add(_ordinates[sourceIndex].Y);
                        yValues.Add(tOrdinates[transformIndex].X);
                        transformIndex += transformStep;
                        sourceIndex += sourceStep;
                    }
                    else if (_ordinates[sourceIndex].X > tOrdinates[transformIndex].Y)
                    {
                        xValues.Add(InterpolateYFromX(tOrdinates[transformIndex].Y, _ordinates[sourceIndex - 1], _ordinates[sourceIndex], xTransform, yTransform)); // GetYfromX(tOrdinates(transformIndex).Y))
                        yValues.Add(tOrdinates[transformIndex].X);
                        transformIndex += transformStep;
                    }
                    else
                    {
                        xValues.Add(_ordinates[sourceIndex].Y);
                        yValues.Add(InterpolateXFromY(_ordinates[sourceIndex].X, tOrdinates[transformIndex - 1], tOrdinates[transformIndex], xTransform, yTransform)); // transform.GetXfromY(_ordinates(sourceIndex).X))
                        sourceIndex += sourceStep;
                    }
                }
            }
            else if (sourceCommonVariableX == false & transformCommonVariableX)
            {
                while (!(sourceIndex < 0 | sourceIndex == _ordinates.Count | transformIndex < 0 | transformIndex == tOrdinates.Count))
                {
                    if (_ordinates[sourceIndex].Y == tOrdinates[transformIndex].X)
                    {
                        xValues.Add(_ordinates[sourceIndex].X);
                        yValues.Add(tOrdinates[transformIndex].Y);
                        transformIndex += transformStep;
                        sourceIndex += sourceStep;
                    }
                    else if (_ordinates[sourceIndex].Y > tOrdinates[transformIndex].X)
                    {
                        xValues.Add(InterpolateXFromY(tOrdinates[transformIndex].X, _ordinates[sourceIndex - 1], _ordinates[sourceIndex], xTransform, yTransform)); // GetXfromY(tOrdinates(transformIndex).X))
                        yValues.Add(tOrdinates[transformIndex].Y);
                        transformIndex += transformStep;
                    }
                    else
                    {
                        xValues.Add(_ordinates[sourceIndex].X);
                        yValues.Add(InterpolateYFromX(_ordinates[sourceIndex].Y, tOrdinates[transformIndex - 1], tOrdinates[transformIndex], xTransform, yTransform)); // .GetYfromX(_ordinates(sourceIndex).Y))
                        sourceIndex += sourceStep;
                    }
                }
            }
            else if (sourceCommonVariableX == false & transformCommonVariableX == false)
            {
                while (!(sourceIndex < 0 | sourceIndex == _ordinates.Count | transformIndex < 0 | transformIndex == tOrdinates.Count))
                {
                    if (_ordinates[sourceIndex].Y == tOrdinates[transformIndex].Y)
                    {
                        xValues.Add(_ordinates[sourceIndex].X);
                        yValues.Add(tOrdinates[transformIndex].X);
                        transformIndex += transformStep;
                        sourceIndex += sourceStep;
                    }
                    else if (_ordinates[sourceIndex].Y > tOrdinates[transformIndex].Y)
                    {
                        xValues.Add(InterpolateXFromY(tOrdinates[transformIndex].Y, _ordinates[sourceIndex - sourceStep], _ordinates[sourceIndex], xTransform, yTransform)); // GetXfromY(tOrdinates(transformIndex).Y))
                        yValues.Add(tOrdinates[transformIndex].X);
                        transformIndex += transformStep;
                    }
                    else
                    {
                        xValues.Add(_ordinates[sourceIndex].X);
                        yValues.Add(InterpolateXFromY(_ordinates[sourceIndex].Y, tOrdinates[transformIndex - transformStep], tOrdinates[transformIndex], xTransform, yTransform)); // transform.GetXfromY(_ordinates(sourceIndex).Y))
                        sourceIndex += sourceStep;
                    }
                }
            }
            // 
            bool strictOnX = StrictX;
            var orderOnX = OrderX;
            if (sourceCommonVariableX == true)
            {
                strictOnX = StrictY;
                orderOnX = OrderY;
            }

            bool strictOnY = transformFunction.StrictX;
            var orderOnY = transformFunction.OrderX;
            if (transformCommonVariableX == true)
            {
                strictOnY = transformFunction.StrictY;
                orderOnY = transformFunction.OrderY;
            }

            return new OrderedPairedData(xValues, yValues, strictOnX, orderOnX, strictOnY, orderOnY);
        }

        /// <summary>
        /// Searches ordinates in the OrderedPairedData collection
        /// for an ordinate based on X-value and returns the zero-based index of the element.
        /// </summary>
        /// <param name="value">The X-value to locate.</param>
        /// <returns>The zero-based index of the ordinate in the collection if item is found; otherwise, a negative number that is the bitwise complement
        /// of the index of the next element that is larger than item or, if there is no larger element, the bitwise complement of the collection count.
        /// </returns>
        public int BinarySearchX(double value)
        {
            int low = 0;
            int high = _ordinates.Count - 1;
            if (OrderX == SortOrder.Ascending)
            {
                while (low <= high)
                {
                    int medianIndex = low + (high - low >> 1);
                    if (_ordinates[medianIndex].X == value)
                        return medianIndex;
                    if (_ordinates[medianIndex].X < value)
                    {
                        low = medianIndex + 1;
                    }
                    else
                    {
                        high = medianIndex - 1;
                    }
                }

                return ~low;
            }
            else
            {
                while (low <= high)
                {
                    int medianIndex = low + (high - low >> 1);
                    if (_ordinates[medianIndex].X == value)
                        return medianIndex;
                    if (_ordinates[medianIndex].X < value)
                    {
                        high = medianIndex - 1;
                    }
                    else
                    {
                        low = medianIndex + 1;
                    }
                }

                return ~low;
            }
        }

        /// <summary>
        /// Searches ordinates in the OrderedPairedData collection
        /// for an ordinate based on Y-value and returns the zero-based index of the element.
        /// </summary>
        /// <param name="value">The Y-value to locate.</param>
        /// <returns>The zero-based index of the ordinate in the collection if item is found; otherwise, a negative number that is the bitwise complement
        /// of the index of the next element that is larger than item or, if there is no larger element, the bitwise complement of the collection count.
        /// </returns>
        public int BinarySearchY(double value)
        {
            int low = 0;
            int high = _ordinates.Count - 1;
            if (OrderY == SortOrder.Ascending)
            {
                while (low <= high)
                {
                    int medianIndex = low + (high - low >> 1);
                    if (_ordinates[medianIndex].Y == value)
                        return medianIndex;
                    if (_ordinates[medianIndex].Y < value)
                    {
                        low = medianIndex + 1;
                    }
                    else
                    {
                        high = medianIndex - 1;
                    }
                }

                return ~low;
            }
            else
            {
                while (low <= high)
                {
                    int medianIndex = low + (high - low >> 1);
                    if (_ordinates[medianIndex].Y == value)
                        return medianIndex;
                    if (_ordinates[medianIndex].Y < value)
                    {
                        high = medianIndex - 1;
                    }
                    else
                    {
                        low = medianIndex + 1;
                    }
                }

                return ~low;
            }
        }

        #endregion

        #region Search Methods

        /// <summary>
        /// Search for the lower bound of the interpolation interval. This method updates whether the values being searched on repeated calls are correlated, 
        /// and saves search value for future use on the next call.  
        /// </summary>
        /// <param name="x">The value to search for.</param>
        public int SearchX(double x)
        {
            int start = 0;
            if (UseSmartSearch)
            {
                start = Xcorrelated ? HuntSearchX(x) : BisectionSearchX(x);
            }
            else
            {
                start = SequentialSearchX(x);
            }
            Xcorrelated = Math.Abs(start - XSearchStart) > XdeltaStart ? false : true;
            XSearchStart = start < 0 || start >= Count ? 0 : start;
            return start;
        }

        /// <summary>
        /// Search for the lower bound of the interpolation interval. This method updates whether the values being searched on repeated calls are correlated, 
        /// and saves search value for future use on the next call.  
        /// </summary>
        /// <param name="y">The value to search for.</param>
        public int SearchY(double y)
        {
            int start = 0;
            if (UseSmartSearch)
            {
                start = Ycorrelated ? HuntSearchY(y) : BisectionSearchY(y);
            }
            else
            {
                start = SequentialSearchY(y);
            }
            Ycorrelated = Math.Abs(start - YSearchStart) > YdeltaStart ? false : true;
            YSearchStart = start < 0 || start >= Count ? 0 : start;
            return start;
        }

        /// <summary>
        /// Searches for the lower bound of the location of a value using a sequential search method.
        /// </summary>
        /// <param name="x">The value to search for.</param>
        public int SequentialSearchX(double x)
        {
            int jl = XSearchStart;
            if ((OrderX == SortOrder.Ascending && x < _ordinates[0].X) ||
                (OrderX == SortOrder.Descending && x > _ordinates[0].X))
            {
                return 0;
            }
            else if ((OrderX == SortOrder.Ascending && x > _ordinates[Count - 1].X) ||
                        (OrderX == SortOrder.Descending && x < _ordinates[Count - 1].X))
            {
                return Count - 2;
            }
            else if ((OrderX == SortOrder.Ascending && x < _ordinates[XSearchStart].X) ||
                        (OrderX == SortOrder.Descending && x > _ordinates[XSearchStart].X))
            {
                jl = 0;
            }
            for (int i = jl; i < Count; i++)
            {
                if ((OrderX == SortOrder.Ascending && x <= _ordinates[i].X) ||
                    (OrderX == SortOrder.Descending && x >= _ordinates[i].X))
                {
                    jl = i - 1;
                    break;
                }
            }
            return jl;
        }

        /// <summary>
        /// Searches for the lower bound of the location of a value using a sequential search method.
        /// </summary>
        /// <param name="y">The value to search for.</param>
        public int SequentialSearchY(double y)
        {
            int jl = YSearchStart;
            if ((OrderY == SortOrder.Ascending && y < _ordinates[0].Y) ||
                (OrderY == SortOrder.Descending && y > _ordinates[0].Y))
            {
                return 0;
            }
            else if ((OrderY == SortOrder.Ascending && y > _ordinates[Count - 1].Y) ||
                        (OrderY == SortOrder.Descending && y < _ordinates[Count - 1].Y))
            {
                return Count - 2;
            }
            else if ((OrderY == SortOrder.Ascending && y < _ordinates[XSearchStart].Y) ||
                        (OrderY == SortOrder.Descending && y > _ordinates[XSearchStart].Y))
            {
                jl = 0;
            }
            for (int i = jl; i < Count; i++)
            {
                if ((OrderY == SortOrder.Ascending && y <= _ordinates[i].Y) ||
                    (OrderY == SortOrder.Descending && y >= _ordinates[i].Y))
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
        public int BisectionSearchX(double x)
        {
            int ju = Count - 1, jm, jl = 0;
            bool ascnd = OrderX == SortOrder.Ascending;
            while (ju - jl > 1)
            {
                jm = (ju + jl) >> 1;
                if (x >= _ordinates[jm].X == ascnd)
                    jl = jm;
                else
                    ju = jm;
            }
            return jl;
        }

        /// <summary>
        /// Searches for the lower bound of the location of a value using a bisection search method.
        /// </summary>
        /// <param name="y">The value to search for.</param>
        public int BisectionSearchY(double y)
        {
            int ju = Count - 1, jm, jl = 0;
            bool ascnd = OrderY == SortOrder.Ascending;
            while (ju - jl > 1)
            {
                jm = (ju + jl) >> 1;
                if (y >= _ordinates[jm].Y == ascnd)
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
        public int HuntSearchX(double x)
        {
            int jl = XSearchStart, jm, ju, inc = 1;
            bool ascnd = OrderX == SortOrder.Ascending;
            if (jl < 0 || jl > Count - 1)
            {
                jl = 0;
                ju = Count - 1;
            }
            else
            {
                if (x >= _ordinates[jl].X == ascnd)
                {
                    for (; ; )
                    {
                        ju = jl + inc;
                        if (ju >= Count - 1) { ju = Count - 1; break; }
                        else if (x < _ordinates[ju].X == ascnd) break;
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
                        else if (x >= _ordinates[jl].X == ascnd) break;
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
                if (x >= _ordinates[jm].X == ascnd)
                    jl = jm;
                else
                    ju = jm;
            }
            return jl;
        }

        /// <summary>
        /// Searches for the lower bound of the location of a value using a hunt search method.
        /// </summary>
        /// <param name="y">The value to search for.</param>
        public int HuntSearchY(double y)
        {
            int jl = YSearchStart, jm, ju, inc = 1;
            bool ascnd = OrderY == SortOrder.Ascending;
            if (jl < 0 || jl > Count - 1)
            {
                jl = 0;
                ju = Count - 1;
            }
            else
            {
                if (y >= _ordinates[jl].Y == ascnd)
                {
                    for (; ; )
                    {
                        ju = jl + inc;
                        if (ju >= Count - 1) { ju = Count - 1; break; }
                        else if (y < _ordinates[ju].Y == ascnd) break;
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
                        else if (y >= _ordinates[jl].Y == ascnd) break;
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
                if (y >= _ordinates[jm].Y == ascnd)
                    jl = jm;
                else
                    ju = jm;
            }
            return jl;
        }

        #endregion

        /// <summary>
        /// Copies the ordinates to an System.Array, starting at a particular System.Array index.
        /// </summary>
        /// <param name="array">The one-dimensional System.Array that is the destination of the elements copied
        /// from. The System.Array must have zero-based indexing.</param>
        /// <param name="arrayIndex">The zero-based index in array at which copying begins.</param>
        public void CopyTo(Ordinate[] array, int arrayIndex)
        {
            _ordinates.CopyTo(array, arrayIndex);
        }

        /// <summary>
        /// Clones the object to a new object.
        /// </summary>
        /// <returns>Clone of target OrderedPairedData collection.</returns>
        public OrderedPairedData Clone()
        {
            return new OrderedPairedData(_ordinates, StrictX, OrderX, StrictY, OrderY);
        }

        /// <summary>
        /// Create an inverted function.
        /// </summary>
        public OrderedPairedData Invert()
        {
            var invertedOrdinates = new Ordinate[_ordinates.Count];
            for (int i = 0; i < _ordinates.Count; i++)
                invertedOrdinates[i] = new Ordinate(_ordinates[i].Y, _ordinates[i].X);
            return new OrderedPairedData(invertedOrdinates, StrictY, OrderY, StrictX, OrderX);
        }



        /// <summary>
        /// Returns an enumerator that iterates through the collection.
        /// </summary>
        /// <returns>An enumerator for the collection.</returns>
        public IEnumerator<Ordinate> GetEnumerator()
        {
            return _ordinates.GetEnumerator();
        }

        /// <summary>
        /// Returns an enumerator that iterates through the collection.
        /// </summary>
        /// <returns>An enumerator for the collection.</returns>
        IEnumerator IEnumerable.GetEnumerator()
        {
            return _ordinates.GetEnumerator();
        }

        /// <summary>
        /// Tests for numerical equality between two OrderedPairedData collection.
        /// </summary>
        /// <param name="left">OrderedPairedData object to the left of the equality operator.</param>
        /// <param name="right">OrderedPairedData object to the right of the equality operator.</param>
        /// <returns>True if two objects are numerically equal; otherwise, False.</returns>
        public static bool operator ==(OrderedPairedData left, OrderedPairedData right)
        {
            // Check for null arguments. Keep in mind null == null
            if (left is null && right is null)
            {
                return true;
            }
            else if (left is null)
            {
                return false;
            }
            else if (right is null)
            {
                return false;
            }
            // I don't think this is possible
            if ((left._ordinates == null) && (right._ordinates == null))
            {
                return true;
            }
            else if (left._ordinates == null)
            {
                return false;
            }
            else if (right._ordinates == null)
            {
                return false;
            }
            if (left.Count != right.Count) return false;
            for (int i = 0; i < left._ordinates.Count; i++)
            {
                if (Math.Abs(left._ordinates[i].X - right._ordinates[i].X) > 0.0000000000001d)
                    return false;
                if (Math.Abs(left._ordinates[i].Y - right._ordinates[i].Y) > 0.0000000000001d)
                    return false;
            }
            return true;
        }

        /// <summary>
        /// Tests for numerical inequality between two OrderedPairedData collection.
        /// </summary>
        /// <param name="left">OrderedPairedData object to the left of the inequality operator.</param>
        /// <param name="right">OrderedPairedData object to the right of the inequality operator.</param>
        /// <returns>True if two objects are not numerically equal; otherwise, False.</returns>
        public static bool operator !=(OrderedPairedData left, OrderedPairedData right)
        {
            return !(left == right);
        }


        #region Line Simplification


        public OrderedPairedData DouglasPeuckerSimplify(double tolerance)
        {
            List<Ordinate> ordinates = new List<Ordinate>();
            var ints = DouglasPeuckerReduction(tolerance);
            for (int i = 0; i < ints.Count; i++)
                ordinates.Add(new Ordinate(_ordinates[ints[i]].X, _ordinates[ints[i]].Y));

            return new OrderedPairedData(ordinates, StrictX, OrderX, StrictY, OrderY);
        }


        /// <summary>
        /// Returns a list of indices to keep for the simplified x and y values of a line.
        /// Source: http://www.codeproject.com/Articles/18936/A-Csharp-Implementation-of-Douglas-Peucker-Line-Ap
        /// </summary>
        /// <param name="xValues">x-values of the line.</param>
        /// <param name="yValues">y-values of the line.</param>
        /// <param name="tolerance">Tolerance to remove points. Higher tolerance will remove more points.</param>
        /// <returns></returns>
        private List<int> DouglasPeuckerReduction(double tolerance)
        {

            // http://www.codeproject.com/Articles/18936/A-Csharp-Implementation-of-Douglas-Peucker-Line-Ap

            List<int> reducedPointIndexes = new List<int>();
            if (_ordinates == null) { return reducedPointIndexes; }
            if (_ordinates.Count() < 3 | tolerance <= 0)
            {
                for (int i = 0; i <= _ordinates.Count() - 1; i++) { reducedPointIndexes.Add(i); }
                return reducedPointIndexes;
            }

            int firstPoint = 0;
            int lastPoint = _ordinates.Count() - 1;
            reducedPointIndexes = new List<int> { firstPoint, lastPoint };
            DouglasPeuckerReduction(firstPoint, lastPoint, tolerance, ref reducedPointIndexes);
            reducedPointIndexes.Sort();
            return reducedPointIndexes;

        }

        private void DouglasPeuckerReduction(int firstPoint, int lastPoint, double tolerance, ref List<int> pointIndexesToKeep)

        {

            double maxDistance = 0;
            int indexFarthest = 0;
            double distance;

            if (firstPoint != (lastPoint - 1))
            {
                for (int index = firstPoint; index <= lastPoint - 1; index++)
                {
                    distance = PerpendicularDistance(_ordinates[firstPoint].X, _ordinates[firstPoint].Y, _ordinates[lastPoint].X, _ordinates[lastPoint].Y, _ordinates[index].X, _ordinates[index].Y);
                    if (distance > maxDistance)
                    {
                        maxDistance = distance;
                        indexFarthest = index;
                    }
                }
            }

            if (maxDistance > tolerance & indexFarthest != 0)
            {
                pointIndexesToKeep.Add(indexFarthest);
                DouglasPeuckerReduction(firstPoint, indexFarthest, tolerance, ref pointIndexesToKeep);
                DouglasPeuckerReduction(indexFarthest, lastPoint, tolerance, ref pointIndexesToKeep);
            }

        }

        private double PerpendicularDistance(double aX, double aY, double bX, double bY, double cX, double cY)
        {
            // Area = |(1/2)(x1y2 + x2y3 + x3y1 - x2y1 - x3y2 - x1y3)|   *Area of triangle
            // Base = v((x1-x2)²+(x1-x2)²)                               *Base of Triangle*
            // Area = .5*Base*H                                          *Solve for height
            // Height = Area/.5/Base

            double area = Math.Abs((aX * bY + bX * cY + cX * aY - bX * aY - cX * bY - aX * cY) * 0.5);
            double triangleBase = Math.Pow(Math.Pow(aX - bX, 2) + Math.Pow(aY - bY, 2), 0.5);
            return area * 2 / triangleBase;
        }

        /// <summary>
        /// Upgrade when we implement the priority queue in .NET 6
        /// </summary>
        /// <param name="numToKeep"></param>
        /// <returns></returns>
        public OrderedPairedData VisvaligamWhyattSimplify(int numToKeep)
        {
            // http://bost.ocks.org/mike/simplify/

            List<Ordinate> ordinates = new List<Ordinate>(_ordinates);

            int removeLimit = ordinates.Count - numToKeep;
            int minIndex;
            double minArea;
            double tmpArea;


            for (int i = 0; i < removeLimit; i++)
            {
                minIndex = 1;

                minArea = TriangleArea(ordinates[0], ordinates[1], ordinates[2]);



                for (int j = 2; j <= ordinates.Count - 2; j++)
                {
                    tmpArea = TriangleArea(ordinates[j - 1], ordinates[j], ordinates[j + 1]);

                    if (tmpArea < minArea)
                    {
                        minIndex = j;
                        minArea = tmpArea;
                    }
                }

                ordinates.RemoveAt(minIndex);
            }

            return new OrderedPairedData(ordinates, StrictX, OrderX, StrictY, OrderY);
        }


        //private bool CompareTriangles(Triangle A, Triangle B)
        //{
        //    return A.Area < B.Area;
        //}
        //private class Triangle
        //{
        //    public int[] Indices = new int[3];
        //    public double Area;
        //    public Triangle Prev;
        //    public Triangle Next;
        //}

        private double TriangleArea(Ordinate point1, Ordinate point2, Ordinate point)
        {
            return Math.Abs((point1.X * point2.Y + point2.X * point.Y + point.X * point1.Y - point2.X * point1.Y - point.X * point2.Y - point1.X * point.Y) * 0.5);
        }

        //public OrderedPairedData VWSimplify(int length)
        //{
        //    List<Ordinate> ordinates = new List<Ordinate>(_ordinates);


        //}

        public OrderedPairedData LangSimplify(double tolerance, int lookAhead)
        {
            if (_ordinates == null | lookAhead <= 1 | tolerance <= 0)
                return this;

            List<Ordinate> ordinates = new List<Ordinate>();

            int count = _ordinates.Count;
            int offset;
            if (lookAhead > count - 1)
                lookAhead = count - 1;
            ordinates.Add(_ordinates[0]);

            for (int i = 0; i < count; i++)
            {
                if (i + lookAhead > count)
                    lookAhead = count - i - 1;

                offset = RecursiveTolerance(i, lookAhead, tolerance);

                if ((offset > 0) & (i + offset < _ordinates.Count))
                {
                    ordinates.Add(new Ordinate(_ordinates[i + offset].X, _ordinates[i + offset].Y));
                    i += offset - 1;
                }
            }

            return new OrderedPairedData(ordinates, StrictX, OrderX, StrictY, OrderY);
        }

        private int RecursiveTolerance(int i, int lookAhead, double tolerance)
        {
            int n = lookAhead;
            double angle;

            var cp = _ordinates[i];

            if (i + n < _ordinates.Count)
            {
                var v1 = new Ordinate(_ordinates[i + n].X - cp.X, _ordinates[i + n].Y - cp.Y);

                for (int j = 1; j <= n; j++)
                {
                    var clp = _ordinates[i + j];

                    var v2 = new Ordinate(clp.X - cp.X, clp.Y - cp.Y);

                    angle = Math.Acos((v1.X * v2.X + v1.Y * v2.Y) / (Math.Sqrt(v1.Y * v1.Y + v1.X * v1.X) * Math.Sqrt(v2.Y * v2.Y + v2.X * v2.X)));

                    if (double.IsNaN(angle) || double.IsInfinity(angle)) angle = 0;

                    double lH = Math.Sqrt((clp.X - cp.X) * (clp.X - cp.X) + (clp.Y - cp.Y) * (clp.Y - cp.Y));

                    if (Math.Sin(angle) * lH >= tolerance)
                    {
                        n -= 1;

                        if (n > 0)
                            return RecursiveTolerance(i, n, tolerance);
                        else
                            return 0;
                    }
                }
            }

            return n;
        }






        #endregion

        /// <summary>
        /// Converts the ordered paired data set to an XElement for saving to xml.
        /// </summary>
        /// <returns>An XElement representation of the data.</returns>
        public XElement SaveToXElement()
        {
            var result = new XElement(nameof(OrderedPairedData));
            // 
            // Order
            result.SetAttributeValue(nameof(StrictX), StrictX.ToString());
            result.SetAttributeValue(nameof(StrictY), StrictY.ToString());
            // Get Strictness
            result.SetAttributeValue(nameof(OrderX), OrderX.ToString());
            result.SetAttributeValue(nameof(OrderY), OrderY.ToString());
            // 
            var curveElement = new XElement("Ordinates");
            for (int i = 0; i < Count; i++) { curveElement.Add(this[i].ToXElement()); }
            // 
            result.Add(curveElement);
            return result;
        }

    }
}