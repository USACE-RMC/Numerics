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
using System.Collections;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Globalization;
using System.Linq;
using System.Xml.Linq;
using Numerics.Distributions;

namespace Numerics.Data
{
    /// <summary>
    /// Class designed to store xy data that is ordered for both the x and y values. Here y is uncertain and represented 
    /// as a continuous distribution.
    /// </summary>
    /// <remarks>
    ///     <b> Authors:</b>
    ///     Woodrow Fields, USACE Risk Management Center, woodrow.l.fields@usace.army.mil
    /// </remarks> 
    [Serializable]
    public class UncertainOrderedPairedData : IList<UncertainOrdinate>, INotifyCollectionChanged
    {
        #region Members
        // 
        private List<UncertainOrdinate> _uncertainOrdinates;
        private bool _strictX;
        private bool _strictY;
        private SortOrder _orderX;
        private SortOrder _orderY;
        private bool _isValid = true;
        private bool _allowDifferentDistributionTypes = false;

        /// <summary>
        /// Gets or sets whether to allow different distribution types in the list. 
        /// </summary>
        public bool AllowDifferentDistributionTypes
        {
            get { return _allowDifferentDistributionTypes; }
            set
            {
                if (_allowDifferentDistributionTypes != value)
                {
                    _allowDifferentDistributionTypes = value;
                    Validate();
                }
            }
        }

        /// <summary>
        /// Determines if the list is valid. 
        /// </summary>
        public bool IsValid
        {
            get { return _isValid; }
        }

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
        /// Gets the univariate distribution type. 
        /// </summary>
        public UnivariateDistributionType Distribution { get; private set; }

        /// <summary>
        /// Gets the count of the number of UncertainOrdinates currently stored
        /// </summary>
        public int Count => _uncertainOrdinates.Count;

        /// <summary>
        /// Boolean that dictates if an object is read-only or can be modified
        /// </summary>
        public bool IsReadOnly => false;

        /// <summary>
        /// Boolean of whether or not to invoke anything on CollectionChanged
        /// </summary>
        public bool SupressCollectionChanged { get; set; } = false;

        /// <summary>
        /// Handles the event of CollectionChanged
        /// </summary>
        public event NotifyCollectionChangedEventHandler CollectionChanged;

        #endregion

        #region Construction

        /// <summary>
        /// Create empty instance of the uncertain ordered paired data class.
        /// </summary>
        /// <param name="strictOnX">Boolean indicating if x values are strictly increasing/decreasing. True means x values cannot be equal.</param>
        /// <param name="xOrder">Order of the x values.</param>
        /// <param name="strictOnY">Boolean indicating if x values are strictly increasing/decreasing. True means x values cannot be equal.</param>
        /// <param name="yOrder">Order of the y values.</param>
        /// <param name="distributionType">The distribution type of the y values</param>
        /// <param name="capacity"> Optional, capacity of the collection.</param>
        public UncertainOrderedPairedData(bool strictOnX, SortOrder xOrder, bool strictOnY, SortOrder yOrder, UnivariateDistributionType distributionType, int capacity = 0)
        {
            if (capacity > 0)
                _uncertainOrdinates = new List<UncertainOrdinate>(capacity);
            else
                _uncertainOrdinates = new List<UncertainOrdinate>();

            Distribution = distributionType;
            _strictX = strictOnX;
            _strictY = strictOnY;
            _orderX = xOrder;
            _orderY = yOrder;
        }

        /// <summary>
        /// Create an instance of the uncertain ordered paired data class with defined uncertain ordinate data.
        /// </summary>
        /// <param name="xData">Uncertain ordinate x values.</param>
        /// <param name="yData">Uncertain ordinate y values.</param>
        /// <param name="strictOnX">Boolean indicating if x values are strictly increasing/decreasing. True means x values cannot be equal.</param>
        /// <param name="xOrder">Order of the x values.</param>
        /// <param name="strictOnY">Boolean indicating if x values are strictly increasing/decreasing. True means x values cannot be equal.</param>
        /// <param name="yOrder">Order of the y values.</param>
        /// <param name="distributionType">The distribution type of the y values</param>
        public UncertainOrderedPairedData(IList<double> xData, IList<UnivariateDistributionBase> yData, bool strictOnX, SortOrder xOrder, bool strictOnY, SortOrder yOrder, UnivariateDistributionType distributionType)
        {
            if (xData.Count != yData.Count)
                throw new Exception("Number of X values and Y values must be the same.");
            Distribution = distributionType;
            _strictX = strictOnX;
            _strictY = strictOnY;
            _orderX = xOrder;
            _orderY = yOrder;
            _uncertainOrdinates = new List<UncertainOrdinate>(xData.Count);
            for (int i = 0; i < xData.Count; i++)
                _uncertainOrdinates.Add(new UncertainOrdinate(xData[i], yData[i].Clone()));
            Validate();
        }

        /// <summary>
        /// Create an instance of the uncertain ordered paired data class with defined uncertain ordinate data.
        /// </summary>
        /// <param name="data">Uncertain ordinate values.</param>
        /// <param name="strictOnX">Boolean indicating if x values are strictly increasing/decreasing. True means x values cannot be equal.</param>
        /// <param name="xOrder">Order of the x values.</param>
        /// <param name="strictOnY">Boolean indicating if x values are strictly increasing/decreasing. True means x values cannot be equal.</param>
        /// <param name="yOrder">Order of the y values.</param>
        /// <param name="distributionType">The distribution type of the y values</param>
        public UncertainOrderedPairedData(IList<UncertainOrdinate> data, bool strictOnX, SortOrder xOrder, bool strictOnY, SortOrder yOrder, UnivariateDistributionType distributionType)
        {
            Distribution = distributionType;
            _strictX = strictOnX;
            _strictY = strictOnY;
            _orderX = xOrder;
            _orderY = yOrder;
            _uncertainOrdinates = new List<UncertainOrdinate>(data.Count);
            for (int i = 0; i < data.Count; i++)
                _uncertainOrdinates.Add(new UncertainOrdinate(data[i].X, data[i].Y.Clone()));
            Validate();
        }

        /// <summary>
        /// Helper function for the Clone() method that creates an instance of the uncertain ordered paired data class with defined uncertain ordinate data.
        /// </summary>
        /// <param name="data">Uncertain ordinate values.</param>
        /// <param name="strictOnX">Boolean indicating if x values are strictly increasing/decreasing. True means x values cannot be equal.</param>
        /// <param name="xOrder">Order of the x values.</param>
        /// <param name="strictOnY">Boolean indicating if x values are strictly increasing/decreasing. True means x values cannot be equal.</param>
        /// <param name="yOrder">Order of the y values.</param>
        /// <param name="distributionType">The distribution type of the y values</param>
        /// <param name="dataValid">Boolean of if the data is valid.</param>
        private UncertainOrderedPairedData(IList<UncertainOrdinate> data, bool strictOnX, SortOrder xOrder, bool strictOnY, SortOrder yOrder, UnivariateDistributionType distributionType, bool dataValid)
        {
            Distribution = distributionType;
            _strictX = strictOnX;
            _strictY = strictOnY;
            _orderX = xOrder;
            _orderY = yOrder;
            _uncertainOrdinates = new List<UncertainOrdinate>(data.Count);
            for (int i = 0; i < data.Count; i++)
                _uncertainOrdinates.Add(new UncertainOrdinate(data[i].X, data[i].Y.Clone()));

            _isValid = dataValid;
        }

        /// <summary>
        /// Create a new instance of the uncertain ordered paired data class from an XElement XML object.
        /// </summary>
        /// <param name="el">The XElement the UncertainOrderPairedData object is being created from.</param>
        public UncertainOrderedPairedData(XElement el)
        {
            // Get Order
            if (el.Attribute("X_Strict") != null)
                bool.TryParse(el.Attribute("X_Strict").Value, out _strictX);
            if (el.Attribute("Y_Strict") != null)
                bool.TryParse(el.Attribute("Y_Strict").Value, out _strictY);
            // Get Strictness
            if (el.Attribute("X_Order") != null)
                Enum.TryParse(el.Attribute("X_Order").Value, out _orderX);
            if (el.Attribute("Y_Order") != null)
                Enum.TryParse(el.Attribute("Y_Order").Value, out _orderY);
            // Distribution type
            Distribution = UnivariateDistributionType.Deterministic;
            if (el.Attribute("Distribution") != null)
            {
                var argresult = Distribution;
                Enum.TryParse(el.Attribute("Distribution").Value, out argresult);
                Distribution = argresult;
            }
            // new prop

            if (el.Attribute(nameof(AllowDifferentDistributionTypes)) != null)
            {
                bool.TryParse(el.Attribute(nameof(AllowDifferentDistributionTypes)).Value, out _allowDifferentDistributionTypes);
                // Get Ordinates
                var curveEl = el.Element("Ordinates");
                _uncertainOrdinates = new List<UncertainOrdinate>();
                foreach (XElement ord in curveEl.Elements(nameof(UncertainOrdinate)))
                    _uncertainOrdinates.Add(new UncertainOrdinate(ord));
            }
            else
            {
                var curveEl = el.Element("Ordinates");
                var xData = new List<double>();
                var yData = new List<UnivariateDistributionBase>();
                if (curveEl != null)
                {
                    foreach (XElement o in curveEl.Elements("Ordinate"))
                    {
                        double.TryParse(o.Attribute("X").Value, NumberStyles.Any, CultureInfo.InvariantCulture, out var xout);
                        xData.Add(xout);
                        var dist = UnivariateDistributionFactory.CreateDistribution(Distribution);
                        var props = dist.GetParameterPropertyNames;
                        var paramVals = new double[(props.Count())];
                        for (int i = 0; i < props.Count(); i++)
                        {
                            double.TryParse(o.Attribute(props[i]).Value, NumberStyles.Any, CultureInfo.InvariantCulture, out var result);
                            paramVals[i] = result;
                        }

                        dist.SetParameters(paramVals);
                        yData.Add(dist);
                    }
                }
                // 
                _uncertainOrdinates = new List<UncertainOrdinate>();
                for (int i = 0; i < xData.Count; i++)
                    _uncertainOrdinates.Add(new UncertainOrdinate(xData[i], yData[i].Clone()));
            }

            Validate();
        }

        #endregion

        #region Methods

        /// <summary>
        /// Resets the event handling of the CollectionChanged event
        /// </summary>
        public void RaiseCollectionChangedReset()
        {
            CollectionChanged?.Invoke(this, new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Reset));
        }

        /// <summary>
        /// This function samples a curve using the probability axis of the continuous distribution in the Y value for each x ordinate, which will result in a new paired data curve.
        /// </summary>
        /// <param name="probability">A value between 0 and 1 representing the value to sample from each y continuous distribution</param>
        /// <returns>A sampled curve provided as an OrderedPairedData object.</returns>
        public OrderedPairedData CurveSample(double probability)
        {
            if (probability < 0d)
                probability = 0d;
            if (probability > 1d)
                probability = 1d;
            var result = new Ordinate[_uncertainOrdinates.Count];
            for (int i = 0; i < _uncertainOrdinates.Count; i++)
                result[i] = _uncertainOrdinates[i].GetOrdinate(probability);
            return new OrderedPairedData(result, StrictX, OrderX, StrictY, OrderY);
        }

        /// <summary>
        /// This function samples a curve using the mean of the Y value for each x ordinate, which will result in a new paired data curve.
        /// </summary>
        /// <returns>A sampled curve provided as an OrderedPairedData object.</returns>
        public OrderedPairedData CurveSample()
        {
            var result = new Ordinate[_uncertainOrdinates.Count];
            for (int i = 0; i < _uncertainOrdinates.Count; i++)
                result[i] = _uncertainOrdinates[i].GetOrdinate();
            return new OrderedPairedData(result, StrictX, OrderX, StrictY, OrderY);
        }

        /// <summary>
        /// Creates a deep clone of the object to a new object.
        /// </summary>
        /// <returns>Deep clone of target UncertainOrderedPairedData collection.</returns>
        public UncertainOrderedPairedData Clone()
        {
            return new UncertainOrderedPairedData(_uncertainOrdinates, StrictX, OrderX, StrictY, OrderY, Distribution, _isValid);
        }

        /// <summary>
        /// Method to determine validity of all ordinates in the collection.
        /// </summary>
        public void Validate()
        {
            if (SupressCollectionChanged == true) return;
            _isValid = true; // innocent until proven guilty
            for (int i = 0; i < _uncertainOrdinates.Count; i++)
            {
                if (OrdinateValid(i, false) == false)
                {
                    _isValid = false;
                    return;
                }
            }
        }

        /// <summary>
        /// Determines the validity of a specific ordinate at a specified index.
        /// </summary>
        /// <param name="index">The index of ordinate item.</param>
        /// <param name="lookBackward">Optional parameter to also look backward when determining validity of the ordinate. This parameter is included as an 
        /// optimization for situations where looking at the previous ordinate is not required.</param>
        /// <returns>A boolean that indicates if the ordinate at the specified index is valid or not within the dataset.</returns>
        private bool OrdinateValid(int index, bool lookBackward = true)
        {
            if (index < 0 || index > _uncertainOrdinates.Count - 1)
                return false;
            // Look backward
            if (lookBackward == true & index > 0)
            {
                if (_uncertainOrdinates[index].OrdinateValid(_uncertainOrdinates[index - 1], StrictX, StrictY, OrderX, OrderY, false, AllowDifferentDistributionTypes) == false)
                    return false;
            }
            // Look forward
            if (index < _uncertainOrdinates.Count - 1)
            {
                if (_uncertainOrdinates[index].OrdinateValid(_uncertainOrdinates[index + 1], StrictX, StrictY, OrderX, OrderY, true, AllowDifferentDistributionTypes) == false)
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
            if (_isValid == true)
                return result;
            for (int i = 0; i < _uncertainOrdinates.Count - 1; i++)
                // Look forward
                result.AddRange(_uncertainOrdinates[i].OrdinateErrors(_uncertainOrdinates[i + 1], StrictX, StrictY, OrderX, OrderY, true, AllowDifferentDistributionTypes));
            result.AddRange(_uncertainOrdinates.Last().OrdinateErrors());
            return result;
        }

        /// <summary>
        /// Tests for numerical equality between two UncertainOrderedPairedData collections.
        /// </summary>
        /// <param name="left">UncertainOrderedPairedData object to the left of the equality operator.</param>
        /// <param name="right">UncertainOrderedPairedData object to the right of the equality operator.</param>
        /// <returns>True if two objects are numerically equal; otherwise, False.</returns>
        public static bool operator ==(UncertainOrderedPairedData left, UncertainOrderedPairedData right)
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
            if ((left._uncertainOrdinates == null) && (right._uncertainOrdinates == null))
            {
                return true;
            }
            else if (left._uncertainOrdinates == null)
            {
                return false;
            }
            else if (right._uncertainOrdinates == null)
            {
                return false;
            }
            if (left.Count != right.Count)
                return false;
            for (int i = 0, loopTo = left.Count - 1; i <= loopTo; i++)
            {
                if (left._uncertainOrdinates[i].X != right._uncertainOrdinates[i].X)
                    return false;
                if (left._uncertainOrdinates[i].Y == right._uncertainOrdinates[i].Y == false)
                    return false;
            }
            return true;
        }

        /// <summary>
        /// Tests for numerical inequality between two UncertainOrderedPairedData collections.
        /// </summary>
        /// <param name="left">UncertainOrderedPairedData object to the left of the inequality operator.</param>
        /// <param name="right">UncertainOrderedPairedData object to the right of the inequality operator.</param>
        /// <returns>True if two objects are not numerically equal; otherwise, False.</returns>
        public static bool operator !=(UncertainOrderedPairedData left, UncertainOrderedPairedData right)
        {
            return !(left == right);
        }


        /// <summary>
        /// Get and sets the element at the specific index.
        /// </summary>
        /// <param name="index">The zero-based index of the element to get or set</param>
        public UncertainOrdinate this[int index]
        {
            get { return _uncertainOrdinates[index]; }
            set
            {
                if (_uncertainOrdinates[index] != value)
                {
                    var oldvalue = _uncertainOrdinates[index];
                    _uncertainOrdinates[index] = value;
                    if (_isValid == true)
                    {
                        if (OrdinateValid(index) == false) { _isValid = false; }
                    }
                    else
                    {
                        if (OrdinateValid(index) == true) { Validate(); }
                    }
                    if (SupressCollectionChanged == false)
                    {
                        CollectionChanged?.Invoke(this, new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Replace, value, oldvalue));
                    }
                }
            }
        }

        /// <summary>
        /// Determines the index of a specific ordinate in the collection.
        /// </summary>
        /// <param name="item">The ordinate to locate in the collection.</param>
        /// <returns>The index of ordinate item if found in the list; otherwise, -1.</returns>
        public int IndexOf(UncertainOrdinate item)
        {
            return _uncertainOrdinates.IndexOf(item);
        }

        /// <summary>
        /// Removes the first occurrence of a target ordinate from the collection.
        /// </summary>
        /// <param name="item">The ordinate to remove from the collection.</param>
        /// <returns>True if item was successfully removed from the collection; otherwise, False.</returns>
        public bool Remove(UncertainOrdinate item)
        {
            int itemIndex = IndexOf(item);
            if (itemIndex == -1)
                return false;
            _uncertainOrdinates.RemoveAt(itemIndex);
            Validate();
            if (SupressCollectionChanged == false)
                CollectionChanged?.Invoke(this, new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Remove, item, itemIndex));
            return true;
        }

        /// <summary>
        /// Removes the ordinate from the collection at the specified index.
        /// </summary>
        /// <param name="index">The zero-based index of the item to remove.</param>
        public void RemoveAt(int index)
        {
            var itemRemoved = _uncertainOrdinates[index];
            _uncertainOrdinates.RemoveAt(index);
            Validate();
            if (SupressCollectionChanged == false)
                CollectionChanged?.Invoke(this, new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Remove, itemRemoved, index));
        }

        /// <summary>
        /// Removes a range of ordinates from the collection.
        /// </summary>
        /// <param name="index">The zero-based starting index of the range of elements to remove.</param>
        /// <param name="count">The number of elements to remove.</param>
        public void RemoveRange(int index, int count)
        {
            var itemsRemoved = new UncertainOrdinate[count];
            int counter = 0;
            for (int i = index; i <= index + count - 1; i++)
            {
                itemsRemoved[counter] = _uncertainOrdinates[i];
                counter += 1;
            }

            _uncertainOrdinates.RemoveRange(index, count);
            Validate();
            if (SupressCollectionChanged == false)
                CollectionChanged?.Invoke(this, new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Remove, itemsRemoved, index));
        }

        /// <summary>
        /// Removes a range of ordinates from the collection at specified positions.
        /// </summary>
        /// <param name="rowIndicesToRemove">The zero-based indices of the items to remove.</param>
        public void RemoveRange(int[] rowIndicesToRemove)
        {
            if ((rowIndicesToRemove == null) || rowIndicesToRemove.Count() == 0)
                return;
            // Need to remove from last to first. Can't assume that the array is sorted. need to clone array to ensure the original index array is maintained.
            int[] sortedIndices = (int[])rowIndicesToRemove.Clone();
            Array.Sort(sortedIndices);
            var removedItems = new List<UncertainOrdinate>();
            bool isContinuous = true;
            for (int i = sortedIndices.Count() - 1; i >= 0; i -= 1)
            {
                if (i < sortedIndices.Count() - 1 && sortedIndices[i] != sortedIndices[i + 1] - 1)
                    isContinuous = false;
                removedItems.Add(_uncertainOrdinates[sortedIndices[i]]);
                _uncertainOrdinates.RemoveAt(sortedIndices[i]);
            }
            Validate();
            if (SupressCollectionChanged == false)
            {
                if (isContinuous)
                {
                    CollectionChanged?.Invoke(this, new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Remove, removedItems, sortedIndices[0]));
                }
                else
                {
                    CollectionChanged?.Invoke(this, new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Remove, removedItems));
                }
            }
        }

        /// <summary>
        /// Adds an ordinate to the collection.
        /// </summary>
        /// <param name="item">The ordinate to add to the collection.</param>
        public void Add(UncertainOrdinate item)
        {
            _uncertainOrdinates.Add(item);
            if (OrdinateValid(_uncertainOrdinates.Count - 1) == false)
                _isValid = false;
            if (SupressCollectionChanged == false)
                CollectionChanged?.Invoke(this, new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Add, item, _uncertainOrdinates.Count - 1));
        }

        /// <summary>
        /// Add a series of ordinates to the collection.
        /// </summary>
        /// <param name="items">Items to be added to the collection.</param>
        public void AddRange(IList<UncertainOrdinate> items)
        {
            if ((items == null) || items.Count == 0)
                return;
            int startIndex = _uncertainOrdinates.Count - 1;
            foreach (var ordinate in items)
            {
                _uncertainOrdinates.Add(ordinate);
                if (OrdinateValid(_uncertainOrdinates.Count - 1) == false)
                    _isValid = false;
            }
            if (SupressCollectionChanged == false)
                CollectionChanged?.Invoke(this, new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Add, items.ToList(), startIndex));
        }

        /// <summary>
        /// Inserts an item to the collection of ordinates.
        /// </summary>
        /// <param name="index">The zero-based index at which item should be inserted.</param>
        /// <param name="item">The ordinate to insert into the collection.</param>
        public void Insert(int index, UncertainOrdinate item)
        {
            _uncertainOrdinates.Insert(index, item);
            // only need to set valid state if it is true..if it is already false then inserting can't make it true.
            if (_isValid == true)
            {
                if (OrdinateValid(index) == false)
                    _isValid = false;
            }

            if (SupressCollectionChanged == false)
                CollectionChanged?.Invoke(this, new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Add, item, index));
        }

        /// <summary>
        /// Insert a range of items into the collection of ordinates.
        /// </summary>
        /// <param name="index">The zero-based index at which the items will be inserted.</param>
        /// <param name="items">The ordinates to insert into the collection.</param>
        public void InsertRange(int index, IList<UncertainOrdinate> items)
        {
            if ((items == null) || items.Count == 0)
                return;
            _uncertainOrdinates.InsertRange(index, items);
            for (int i = index; i <= index + items.Count - 1; i++)
            {
                if (_isValid == true)
                {
                    if (OrdinateValid(index) == false)
                        _isValid = false;
                }
            }
            if (SupressCollectionChanged == false)
                CollectionChanged?.Invoke(this, new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Add, items.ToList(), index));
        }

        /// <summary>
        /// Removes all ordinates form the collection and sets count to zero.
        /// </summary>
        public void Clear()
        {
            _uncertainOrdinates.Clear();
            if (SupressCollectionChanged == false)
                CollectionChanged?.Invoke(this, new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Reset));
            _isValid = true;
        }

        /// <summary>
        /// Determines whether the collection contains a specific uncertain ordinate.
        /// </summary>
        /// <param name="item">The uncertain ordinate to locate in the collection.</param>
        /// <returns>True if ordinate is found in the collection; otherwise, false.</returns>
        public bool Contains(UncertainOrdinate item)
        {
            return _uncertainOrdinates.Contains(item);
        }

        /// <summary>
        /// Copies the ordinates to an System.Array, starting at a particular System.Array index.
        /// </summary>
        /// <param name="array">The one-dimensional System.Array that is the destination of the elements copied
        /// from. The System.Array must have zero-based indexing.</param>
        /// <param name="arrayIndex">The zero-based index in array at which copying begins.</param>
        public void CopyTo(UncertainOrdinate[] array, int arrayIndex)
        {
            _uncertainOrdinates.CopyTo(array, arrayIndex);
        }

        /// <summary>
        /// Returns an enumerator that iterates through the collection.
        /// </summary>
        /// <returns>An enumerator for the collection.</returns>
        public IEnumerator<UncertainOrdinate> GetEnumerator()
        {
            return _uncertainOrdinates.GetEnumerator();
        }

        /// <summary>
        /// Returns an enumerator that iterates through the collection.
        /// </summary>
        /// <returns>An enumerator for the collection.</returns>
        IEnumerator IEnumerable.GetEnumerator()
        {
            return _uncertainOrdinates.GetEnumerator();
        }

        /// <summary>
        /// Converts the uncertain ordered paired data set to an XElement for saving to xml.
        /// </summary>
        /// <returns>An XElement representation of the data.</returns>
        public XElement SaveToXElement()
        {
            var result = new XElement("UncertainOrderedPairedData");
            // Order
            result.SetAttributeValue("X_Strict", StrictX.ToString());
            result.SetAttributeValue("Y_Strict", StrictY.ToString());
            // Get Strictness
            result.SetAttributeValue("X_Order", OrderX.ToString());
            result.SetAttributeValue("Y_Order", OrderY.ToString());
            // Distribution type
            result.SetAttributeValue("Distribution", Distribution.ToString());
            result.SetAttributeValue(nameof(AllowDifferentDistributionTypes), AllowDifferentDistributionTypes.ToString());
            // 
            var curveElement = new XElement("Ordinates");
            for (int i = 0; i < Count; i++) { curveElement.Add(this[i].ToXElement()); }
            // 
            result.Add(curveElement);
            return result;
        }

        #endregion
    }
}   