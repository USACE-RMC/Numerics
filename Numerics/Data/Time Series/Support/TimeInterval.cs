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