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
using System.ComponentModel;

namespace Numerics.Data
{

    /// <summary>
    /// Enumeration of available time-series time intervals.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>

    [Serializable]
    public enum TimeInterval
    {
        /// <summary>
        /// A 1 minute time interval.
        /// </summary>
        [Description("1-Min")]
        OneMinute,

        /// <summary>
        /// A 5 minute time interval.
        /// </summary>
        [Description("5-Min")]
        FiveMinute,

        /// <summary>
        /// A 15 minute time interval.
        /// </summary>
        [Description("15-Min")]
        FifteenMinute,

        /// <summary>
        /// A 30 minute time interval.
        /// </summary>
        [Description("30-Min")]
        ThirtyMinute,

        /// <summary>
        /// A 1 hour time interval.
        /// </summary>
        [Description("1-Hr")]
        OneHour,

        /// <summary>
        /// A 6 hour time interval.
        /// </summary>
        [Description("6-Hr")]
        SixHour,

        /// <summary>
        /// A 12 hour time interval.
        /// </summary>
        [Description("12-Hr")]
        TwelveHour,

        /// <summary>
        /// A 1 day time interval.
        /// </summary>
        [Description("1-Day")]
        OneDay,

        /// <summary>
        /// A 7 day, or 1 week time interval.
        /// </summary>
        [Description("7-Day")]
        SevenDay,

        /// <summary>
        /// A 1 month time interval.
        /// </summary>
        [Description("1-Month")]
        OneMonth,

        /// <summary>
        /// A 1 quarter, or 3 month, time interval.
        /// </summary>
        [Description("1-Quarter")]
        OneQuarter,

        /// <summary>
        /// A 1 year time interval.
        /// </summary>
        [Description("1-Year")]
        OneYear,

        /// <summary>
        /// An irregular time interval.
        /// </summary>
        [Description("Irregular")]
        Irregular
    }
}