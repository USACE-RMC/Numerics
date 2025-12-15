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
using System.Collections.Specialized;
using System.ComponentModel;
using System.Globalization;
using System.Linq;

namespace Numerics.Data
{
    /// <summary>
    /// Represents a collection of exceedance probability ordinates used for
    /// plotting and tabulating univariate and multivariate distributions.
    /// </summary>
    /// <remarks>
    /// <para>
    /// This class derives from <see cref="List{Double}"/> so that it can be used
    /// with existing APIs that expect a list. It also implements
    /// <see cref="INotifyCollectionChanged"/> and <see cref="INotifyPropertyChanged"/>
    /// to provide behavior similar to <see cref="System.Collections.ObjectModel.ObservableCollection{T}"/>.
    /// </para>
    /// <para>
    /// The class provides:
    /// </para>
    /// <list type="bullet">
    /// <item><description>Standard default probability ordinates on construction.</description></item>
    /// <item><description>Methods to serialize to and from a delimited string.</description></item>
    /// <item><description>Collection change notifications on key list operations.</description></item>
    /// </list>
    /// <para>
    /// Note: Because <see cref="List{T}"/> methods are not virtual, this class
    /// shadows mutating methods using the <c>new</c> keyword to raise events.
    /// Collection change notifications are guaranteed when the collection is
    /// used as <see cref="ProbabilityOrdinates"/> and the shadowing methods
    /// are called. If the instance is referenced as a plain <see cref="List{Double}"/>,
    /// changes will not raise notifications.
    /// </para>
    /// </remarks>
    public class ProbabilityOrdinates : List<double>, INotifyCollectionChanged, INotifyPropertyChanged
    {
        /// <summary>
        /// The default delimiter used when converting to and from strings.
        /// </summary>
        public const string DefaultDelimiter = "|";

        /// <summary>
        /// Occurs when the collection changes, for example when items are added or removed.
        /// </summary>
        public event NotifyCollectionChangedEventHandler CollectionChanged;

        /// <summary>
        /// Occurs when a property value changes, such as <see cref="Count"/>
        /// or the indexer <c>Item[]</c>.
        /// </summary>
        public event PropertyChangedEventHandler PropertyChanged;

        /// <summary>
        /// Initializes a new instance of the <see cref="ProbabilityOrdinates"/> class
        /// with the standard default probability ordinates.
        /// </summary>
        /// <remarks>
        /// <para>
        /// The default values are:
        /// </para>
        /// <para>
        /// 0.000001, 0.000002, 0.000005, 0.00001, 0.00002, 0.00005,
        /// 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005,
        /// 0.01, 0.02, 0.05, 0.1, 0.2, 0.3,
        /// 0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99.
        /// </para>
        /// </remarks>
        public ProbabilityOrdinates()
        {
            AddDefaults();
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="ProbabilityOrdinates"/> class
        /// populated from the specified sequence of probability values.
        /// </summary>
        /// <param name="probabilities">
        /// Sequence of exceedance probabilities to populate the collection with.
        /// </param>
        /// <exception cref="ArgumentNullException">
        /// Thrown if <paramref name="probabilities"/> is <c>null</c>.
        /// </exception>
        public ProbabilityOrdinates(IEnumerable<double> probabilities)
        {
            if (probabilities == null)
                throw new ArgumentNullException(nameof(probabilities));

            AddRange(probabilities);
        }

        /// <summary>
        /// Returns a string representation of the probability ordinates using
        /// the default delimiter.
        /// </summary>
        /// <returns>
        /// A single string containing all ordinates separated by
        /// <see cref="DefaultDelimiter"/>.
        /// </returns>
        /// <remarks>
        /// <para>
        /// Each value is formatted using the "G17" format string and
        /// <see cref="CultureInfo.InvariantCulture"/>.
        /// </para>
        /// </remarks>
        public override string ToString()
        {
            return ToDelimitedString(DefaultDelimiter);
        }

        /// <summary>
        /// Converts the probability ordinates to a single delimited string.
        /// </summary>
        /// <param name="delimiter">
        /// The delimiter string to insert between values. If <c>null</c> or empty,
        /// the <see cref="DefaultDelimiter"/> is used.
        /// </param>
        /// <returns>
        /// A delimited string containing all ordinates.
        /// </returns>
        public string ToDelimitedString(string delimiter)
        {
            var sep = string.IsNullOrEmpty(delimiter) ? DefaultDelimiter : delimiter;

            return string.Join(
                sep,
                this.Select(p => p.ToString("G17", CultureInfo.InvariantCulture))
            );
        }

        /// <summary>
        /// Clears the collection and replaces its contents with the default
        /// probability ordinates.
        /// </summary>
        public void ResetToDefaults()
        {
            Clear();
            AddDefaults();
        }

        /// <summary>
        /// Clears the collection and populates it from a delimited string of values.
        /// </summary>
        /// <param name="value">
        /// Delimited string containing one or more probability values.
        /// </param>
        /// <param name="delimiter">
        /// The delimiter string separating values in <paramref name="value"/>.
        /// If <c>null</c> or empty, <see cref="DefaultDelimiter"/> is used.
        /// </param>
        /// <remarks>
        /// <para>
        /// This method clears the existing contents and then attempts to parse
        /// each token as a <see cref="double"/> using
        /// <see cref="CultureInfo.InvariantCulture"/>. Empty tokens and tokens that
        /// cannot be parsed are ignored.
        /// </para>
        /// </remarks>
        public void FromDelimitedString(string value, string delimiter = DefaultDelimiter)
        {
            // Use the public Clear method so events are raised.
            Clear();

            if (string.IsNullOrWhiteSpace(value))
                return;

            var sep = string.IsNullOrEmpty(delimiter) ? DefaultDelimiter : delimiter;
            var tokens = value.Split(new[] { sep }, StringSplitOptions.None);

            // Bulk add, then raise Reset for simplicity.
            var toAdd = new List<double>();

            foreach (var token in tokens)
            {
                if (string.IsNullOrWhiteSpace(token))
                    continue;

                if (double.TryParse(
                        token,
                        NumberStyles.Float | NumberStyles.AllowThousands,
                        CultureInfo.InvariantCulture,
                        out var p))
                {
                    toAdd.Add(p);
                }
            }

            if (toAdd.Count > 0)
                AddRange(toAdd);
        }

        /// <summary>
        /// Parses a delimited string of probability ordinates into a new
        /// <see cref="ProbabilityOrdinates"/> instance.
        /// </summary>
        /// <param name="value">
        /// Delimited string containing one or more probability values.
        /// </param>
        /// <param name="delimiter">
        /// The delimiter string separating values in <paramref name="value"/>.
        /// If <c>null</c> or empty, <see cref="DefaultDelimiter"/> is used.
        /// </param>
        /// <returns>
        /// A new <see cref="ProbabilityOrdinates"/> instance populated with any
        /// successfully parsed probability values.
        /// </returns>
        public static ProbabilityOrdinates Parse(string value, string delimiter = DefaultDelimiter)
        {
            var ordinates = new ProbabilityOrdinates();
            ordinates.Clear(); // clears defaults and raises events
            ordinates.FromDelimitedString(value, delimiter);
            return ordinates;
        }

        /// <summary>
        /// Adds the standard default probability ordinates to the collection.
        /// </summary>
        /// <remarks>
        /// <para>
        /// This method does not clear existing values. It appends the default set
        /// to the current contents. Use <see cref="ResetToDefaults"/> to clear
        /// and then add defaults.
        /// </para>
        /// </remarks>
        private void AddDefaults()
        {
            // Use AddRange so events are raised once for the batch.
            AddRange(new[]
            {
                0.000001, 0.000002, 0.000005,
                0.00001,  0.00002,  0.00005,
                0.0001,   0.0002,   0.0005,
                0.001,    0.002,    0.005,
                0.01,     0.02,     0.05,
                0.1,      0.2,      0.3,
                0.5,      0.7,      0.8,
                0.9,      0.95,     0.98,
                0.99
            });
        }

        /// <summary>
        /// Validates the current state of the object and reports any issues found.
        /// </summary>
        /// <returns>
        /// A tuple containing:
        /// <list type="bullet">
        /// <item>
        /// <description><c>IsValid</c>: <c>true</c> if the object passes all validation checks; otherwise <c>false</c>.</description>
        /// </item>
        /// <item>
        /// <description><c>ValidationMessages</c>: a list of messages describing any validation errors or warnings. This list is empty when the object is valid.</description>
        /// </item>
        /// </list>
        /// </returns>
        public (bool IsValid, List<string> ValidationMessages) Validate()
        {
            bool isValid = true;
            var messages = new List<string>();

            // No ordinates at all
            if (Count == 0)
            {
                isValid = false;
                messages.Add("At least one exceedance probability must be specified.");
                return (isValid, messages);
            }

            for (int i = 0; i < Count; i++)
            {
                double p = this[i];

                // Range check
                if (p <= 0.0 || p >= 1.0)
                {
                    isValid = false;
                    messages.Add(
                        $"Exceedance probability at index {i} has value {p.ToString("G17", CultureInfo.InvariantCulture)}, " +
                        "but it must be strictly between 0 and 1.");
                    break;
                }

                // Strictly increasing check
                if (i > 0 && p <= this[i - 1])
                {
                    double prev = this[i - 1];
                    messages.Add(
                        "Exceedance probabilities must be strictly increasing. " +
                        $"Value at index {i - 1} is {prev.ToString("G17", CultureInfo.InvariantCulture)} and " +
                        $"value at index {i} is {p.ToString("G17", CultureInfo.InvariantCulture)}.");
                    isValid = false;
                    break;
                }
            }

            return (isValid, messages);
        }

        #region Mutating methods with change notifications

        /// <summary>
        /// Adds an item to the end of the collection and raises the appropriate
        /// collection and property change events.
        /// </summary>
        /// <param name="item">The probability ordinate to add.</param>
        public new void Add(double item)
        {
            base.Add(item);
            OnPropertyChanged(nameof(Count));
            OnPropertyChanged("Item[]");
            OnCollectionChanged(
                new NotifyCollectionChangedEventArgs(
                    NotifyCollectionChangedAction.Add,
                    item,
                    Count - 1));
        }

        /// <summary>
        /// Adds the elements of the specified collection to the end of the list
        /// and raises a collection reset event.
        /// </summary>
        /// <param name="collection">The collection whose elements should be added.</param>
        public new void AddRange(IEnumerable<double> collection)
        {
            if (collection == null)
                throw new ArgumentNullException(nameof(collection));

            if (!collection.Any())
                return;

            int oldCount = Count;
            base.AddRange(collection);

            OnPropertyChanged(nameof(Count));
            OnPropertyChanged("Item[]");
            // For simplicity, treat AddRange as a Reset notification.
            OnCollectionChanged(
                new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Reset));
        }

        /// <summary>
        /// Inserts an item into the collection at the specified index and raises
        /// the appropriate change events.
        /// </summary>
        /// <param name="index">
        /// The zero-based index at which the item should be inserted.
        /// </param>
        /// <param name="item">The item to insert.</param>
        public new void Insert(int index, double item)
        {
            base.Insert(index, item);
            OnPropertyChanged(nameof(Count));
            OnPropertyChanged("Item[]");
            OnCollectionChanged(
                new NotifyCollectionChangedEventArgs(
                    NotifyCollectionChangedAction.Add,
                    item,
                    index));
        }

        /// <summary>
        /// Inserts the elements of a collection into the list at the specified index
        /// and raises a collection reset event.
        /// </summary>
        /// <param name="index">
        /// The zero-based index at which the new elements should be inserted.
        /// </param>
        /// <param name="collection">The collection whose elements should be inserted.</param>
        public new void InsertRange(int index, IEnumerable<double> collection)
        {
            if (collection == null)
                throw new ArgumentNullException(nameof(collection));

            if (!collection.Any())
                return;

            base.InsertRange(index, collection);

            OnPropertyChanged(nameof(Count));
            OnPropertyChanged("Item[]");
            OnCollectionChanged(
                new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Reset));
        }

        /// <summary>
        /// Removes the first occurrence of a specific object from the collection
        /// and raises the appropriate change events if an item is removed.
        /// </summary>
        /// <param name="item">The item to remove.</param>
        /// <returns>
        /// <c>true</c> if the item was successfully removed; otherwise <c>false</c>.
        /// </returns>
        public new bool Remove(double item)
        {
            int index = IndexOf(item);
            if (index < 0)
                return false;

            bool removed = base.Remove(item);
            if (removed)
            {
                OnPropertyChanged(nameof(Count));
                OnPropertyChanged("Item[]");
                OnCollectionChanged(
                    new NotifyCollectionChangedEventArgs(
                        NotifyCollectionChangedAction.Remove,
                        item,
                        index));
            }
            return removed;
        }

        /// <summary>
        /// Removes the element at the specified index and raises the appropriate
        /// change events.
        /// </summary>
        /// <param name="index">
        /// The zero-based index of the element to remove.
        /// </param>
        public new void RemoveAt(int index)
        {
            if (index < 0 || index >= Count)
                throw new ArgumentOutOfRangeException(nameof(index));

            double removedItem = this[index];
            base.RemoveAt(index);

            OnPropertyChanged(nameof(Count));
            OnPropertyChanged("Item[]");
            OnCollectionChanged(
                new NotifyCollectionChangedEventArgs(
                    NotifyCollectionChangedAction.Remove,
                    removedItem,
                    index));
        }

        /// <summary>
        /// Removes a range of elements from the list and raises a collection reset event.
        /// </summary>
        /// <param name="index">The zero-based starting index of the range.</param>
        /// <param name="count">The number of elements to remove.</param>
        public new void RemoveRange(int index, int count)
        {
            if (count == 0)
                return;

            base.RemoveRange(index, count);

            OnPropertyChanged(nameof(Count));
            OnPropertyChanged("Item[]");
            OnCollectionChanged(
                new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Reset));
        }

        /// <summary>
        /// Removes all elements from the list and raises the appropriate change events.
        /// </summary>
        public new void Clear()
        {
            if (Count == 0)
                return;

            base.Clear();
            OnPropertyChanged(nameof(Count));
            OnPropertyChanged("Item[]");
            OnCollectionChanged(
                new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Reset));
        }

        /// <summary>
        /// Gets or sets the element at the specified index and raises the appropriate
        /// change events when the value is set.
        /// </summary>
        /// <param name="index">The zero-based index of the element.</param>
        /// <returns>The element at the specified index.</returns>
        public new double this[int index]
        {
            get => base[index];
            set
            {
                if (index < 0 || index >= Count)
                    throw new ArgumentOutOfRangeException(nameof(index));

                double oldItem = base[index];
                base[index] = value;

                OnPropertyChanged("Item[]");
                OnCollectionChanged(
                    new NotifyCollectionChangedEventArgs(
                        NotifyCollectionChangedAction.Replace,
                        newItem: value,
                        oldItem: oldItem,
                        index: index));
            }
        }

        #endregion

        #region Helpers for raising events

        /// <summary>
        /// Raises the <see cref="PropertyChanged"/> event for the specified property name.
        /// </summary>
        /// <param name="propertyName">Name of the property that changed.</param>
        protected virtual void OnPropertyChanged(string propertyName)
        {
            var handler = PropertyChanged;
            if (handler != null)
                handler(this, new PropertyChangedEventArgs(propertyName));
        }

        /// <summary>
        /// Raises the <see cref="CollectionChanged"/> event with the specified arguments.
        /// </summary>
        /// <param name="e">
        /// <see cref="NotifyCollectionChangedEventArgs"/> that describes the change.
        /// </param>
        protected virtual void OnCollectionChanged(NotifyCollectionChangedEventArgs e)
        {
            var handler = CollectionChanged;
            if (handler != null)
                handler(this, e);
        }

        #endregion
    }
}

