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


using System.Collections;
using System.Collections.Specialized;

namespace Numerics.Data
{

    /// <summary>
    /// An abstract series class, which is a collection of series ordinates.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public abstract class Series<TIndex, TValue> : IList<SeriesOrdinate<TIndex, TValue>>, IList, INotifyCollectionChanged
    {
        /// <summary>
        /// Internal list.
        /// </summary>
        protected List<SeriesOrdinate<TIndex, TValue>> _seriesOrdinates = new List<SeriesOrdinate<TIndex, TValue>>();

        /// <inheritdoc/>
        public event NotifyCollectionChangedEventHandler CollectionChanged;

        /// <inheritdoc/>
        public SeriesOrdinate<TIndex, TValue> this[int index]
        {
            get { return _seriesOrdinates[index]; }
            set
            {
                if (_seriesOrdinates[index] != value)
                {
                    var oldvalue = _seriesOrdinates[index];
                    _seriesOrdinates[index] = value;
                    if (SuppressCollectionChanged == false)
                        CollectionChanged?.Invoke(this, new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Replace, value, oldvalue));
                }
            }
        }

        /// <inheritdoc/>
        object IList.this[int index]
        {
            get { return _seriesOrdinates[index]; }
            set
            {
                if (value.GetType() != typeof(SeriesOrdinate<TIndex, TValue>))
                {
                    if (_seriesOrdinates[index] != (SeriesOrdinate<TIndex, TValue>)value)
                    {
                        var oldvalue = _seriesOrdinates[index];
                        _seriesOrdinates[index] = (SeriesOrdinate<TIndex, TValue>)value;
                        if (SuppressCollectionChanged == false)
                            CollectionChanged?.Invoke(this, new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Replace, value, oldvalue));
                    }
                }

            }
        }

        /// <summary>
        /// Suppress collection changed events from firing.
        /// </summary>
        public bool SuppressCollectionChanged { get; set; } = false;

        /// <inheritdoc/>
        public int Count => _seriesOrdinates.Count;

        /// <inheritdoc/>
        public virtual bool IsReadOnly => false;

        /// <inheritdoc/>
        public bool IsFixedSize => false;

        /// <inheritdoc/>
        public object SyncRoot => _seriesOrdinates.Count > 0 ? _seriesOrdinates[0] : null;

        /// <inheritdoc/>
        public bool IsSynchronized => false;

        /// <inheritdoc/>
        public virtual void Add(SeriesOrdinate<TIndex, TValue> item)
        {
            if (item == null) throw new ArgumentNullException(nameof(item));
            _seriesOrdinates.Add(item);
            RaiseCollectionChanged(new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Add, item, _seriesOrdinates.Count - 1));
        }

        /// <inheritdoc/>
        public int Add(object item)
        {
            if (item.GetType() != typeof(SeriesOrdinate<TIndex, TValue>)) { return -1; }
            Add((SeriesOrdinate<TIndex, TValue>)item);
            return _seriesOrdinates.Count - 1;
        }

        /// <inheritdoc/>
        public virtual void Insert(int index, SeriesOrdinate<TIndex, TValue> item)
        {
            if (item == null) throw new ArgumentNullException(nameof(item));
            _seriesOrdinates.Insert(index, item);
            RaiseCollectionChanged(new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Add, item, index));
        }

        /// <inheritdoc/>
        public void Insert(int index, object item)
        {
            if (item.GetType() == typeof(SeriesOrdinate<TIndex, TValue>))
            {
                Insert(index, (SeriesOrdinate<TIndex, TValue>)item);
            }
        }

        /// <inheritdoc/>
        public virtual bool Remove(SeriesOrdinate<TIndex, TValue> item)
        {
            var index = IndexOf(item);
            if (_seriesOrdinates.Remove(item) == true)
            {
                RaiseCollectionChanged(new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Remove, item, index));
                return true;
            }
            return false;
        }

        /// <inheritdoc/>
        public void Remove(object item)
        {
            if (item.GetType() == typeof(SeriesOrdinate<TIndex, TValue>))
            {
                Remove((SeriesOrdinate<TIndex, TValue>)item);
            }
        }

        /// <inheritdoc/>
        public virtual void RemoveAt(int index)
        {
            Remove(_seriesOrdinates[index]);
        }

        /// <summary>
        /// Remove all elements from the collection.
        /// </summary>
        public void Clear()
        {
            bool suppress = SuppressCollectionChanged;
            SuppressCollectionChanged = true;
            for (int i = _seriesOrdinates.Count - 1; i >= 0; i--)
                Remove(_seriesOrdinates[i]);
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
            SuppressCollectionChanged = suppress;
        }

        /// <inheritdoc/>
        public bool Contains(SeriesOrdinate<TIndex, TValue> item)
        {
            return _seriesOrdinates.Contains(item);
        }

        /// <inheritdoc/>
        public bool Contains(object item)
        {
            if (item.GetType() == typeof(SeriesOrdinate<TIndex, TValue>))
            {
                return Contains((SeriesOrdinate<TIndex, TValue>)item);
            }
            return false;
        }

        /// <inheritdoc/>
        public void CopyTo(SeriesOrdinate<TIndex, TValue>[] array, int index)
        {
            _seriesOrdinates.CopyTo(array, index);
        }

        /// <inheritdoc/>
        public void CopyTo(Array array, int index)
        {
            if (array.GetType() == typeof(SeriesOrdinate<TIndex, TValue>[]))
            {
                CopyTo((SeriesOrdinate<TIndex, TValue>[])array, index);
            }
        }

        /// <inheritdoc/>
        public int IndexOf(SeriesOrdinate<TIndex, TValue> item)
        {
            return _seriesOrdinates.IndexOf(item);
        }

        /// <inheritdoc/>
        public int IndexOf(object item)
        {
            if (item.GetType() == typeof(SeriesOrdinate<TIndex, TValue>))
            {
                return _seriesOrdinates.IndexOf((SeriesOrdinate<TIndex, TValue>)item);
            }
            return -1;
        }

        /// <inheritdoc/>
        public IEnumerator<SeriesOrdinate<TIndex, TValue>> GetEnumerator()
        {
            return _seriesOrdinates.GetEnumerator();
        }

        /// <inheritdoc/>
        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }

        /// <summary>
        /// Raise the collection changed event.
        /// </summary>
        /// <param name="args">Change event args.</param>
        protected void RaiseCollectionChanged(NotifyCollectionChangedEventArgs args)
        {
            if (!SuppressCollectionChanged)
                CollectionChanged?.Invoke(this, args);
        }

        /// <summary>
        /// Raise the collection changed reset event.
        /// </summary>
        public void RaiseCollectionChangedReset()
        {
            CollectionChanged?.Invoke(this, new NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction.Reset));
        }

        /// <summary>
        /// Returns the list of series values as a list.
        /// </summary>
        public List<TValue> ValuesToList()
        {
            return _seriesOrdinates.Select(x => x.Value).ToList();
        }

        /// <summary>
        /// Returns the list of series values as an array.
        /// </summary>
        public TValue[] ValuesToArray()
        {
            return _seriesOrdinates.Select(x => x.Value).ToArray();
        }

        /// <summary>
        /// Returns the list of series indexes as a list.
        /// </summary>
        public List<TIndex> IndexesToList()
        {
            return _seriesOrdinates.Select(x => x.Index).ToList();
        }

        /// <summary>
        /// Returns the list of series indexes as an array.
        /// </summary>
        public TIndex[] IndexesToArray()
        {
            return _seriesOrdinates.Select(x => x.Index).ToArray();
        }

    }
}

