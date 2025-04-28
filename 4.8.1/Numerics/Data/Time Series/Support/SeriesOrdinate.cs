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
using System.ComponentModel;

namespace Numerics.Data
{

    /// <summary>
    /// A series ordinate.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class SeriesOrdinate<TIndex, TValue> : INotifyPropertyChanged, IEquatable<SeriesOrdinate<TIndex, TValue>>
    {
        /// <summary>
        /// Constructs a new series ordinate.
        /// </summary>
        public SeriesOrdinate() { }

        /// <summary>
        /// Constructs a new series ordinate. 
        /// </summary>
        /// <param name="index">The ordinate index.</param>
        /// <param name="value">The ordinate value.</param>
        public SeriesOrdinate(TIndex index, TValue value)
        {
            _index = index;
            _value = value;
        }

        /// <summary>
        /// Protected index property.
        /// </summary>
        protected TIndex _index;

        /// <summary>
        /// Protected value property.
        /// </summary>
        protected TValue _value;

        /// <inheritdoc/>
        public event PropertyChangedEventHandler PropertyChanged;

        /// <summary>
        /// The index of the series ordinate.
        /// </summary>
        public virtual TIndex Index
        {
            get { return _index; }
            set
            {
                if (!EqualityComparer<TIndex>.Default.Equals(_index, value))
                {
                    _index = value;
                    RaisePropertyChanged(nameof(Index));
                }
            }
        }

        /// <summary>
        /// The value of the time-series ordinate.
        /// </summary>
        public virtual TValue Value
        {
            get { return _value; }
            set
            {
                if (!EqualityComparer<TValue>.Default.Equals(_value, value))
                {
                    _value = value;
                    RaisePropertyChanged(nameof(Value));
                }
            }
        }

        /// <inheritdoc/>
        public bool Equals(SeriesOrdinate<TIndex, TValue> other)
        {
            if (ReferenceEquals(other, null)) return false;
            if (ReferenceEquals(this, other)) return true;
            return EqualityComparer<TIndex>.Default.Equals(_index, other._index)
                && EqualityComparer<TValue>.Default.Equals(_value, other._value);
        }

        /// <inheritdoc/>
        public override bool Equals(object obj) => Equals(obj as SeriesOrdinate<TIndex, TValue>);

        /// <summary>
        /// Equality operator overload. 
        /// </summary>
        /// <param name="left">The first SeriesOrdinate object to compare.</param>
        /// <param name="right">The second SeriesOrdinate object to compare/</param>
        /// <returns>True of the two SeriesOrdinate objects are equal and false otherwise.</returns>
        public static bool operator ==(SeriesOrdinate<TIndex, TValue> left, SeriesOrdinate<TIndex, TValue> right) => EqualityComparer<SeriesOrdinate<TIndex, TValue>>.Default.Equals(left, right);


        /// <summary>
        /// Inequality operator overload. 
        /// </summary>
        /// <param name="left">The first SeriesOrdinate object to compare.</param>
        /// <param name="right">The second SeriesOrdinate object to compare/</param>
        /// <returns>True of the two SeriesOrdinate objects are NOT equal and false otherwise.</returns>
        public static bool operator !=(SeriesOrdinate<TIndex, TValue> left, SeriesOrdinate<TIndex, TValue> right) => !(left == right);
        
        /// <inheritdoc/>
        public override int GetHashCode()
        {
            unchecked
            {
                int hash = 17;
                hash = hash * 23 + EqualityComparer<TIndex>.Default.GetHashCode(_index);
                hash = hash * 23 + EqualityComparer<TValue>.Default.GetHashCode(_value);
                return hash;
            }
        }

        /// <summary>
        /// Raise property changed event.
        /// </summary>
        /// <param name="propertyName">Name of property that changed.</param>
        protected void RaisePropertyChanged(string propertyName)
        {
            PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(propertyName));
        }

        /// <summary>
        /// Returns a copy of the series ordinate.
        /// </summary>
        public SeriesOrdinate<TIndex, TValue> Clone()
        {
            return new SeriesOrdinate<TIndex, TValue>(Index, Value);
        }

    }
}
