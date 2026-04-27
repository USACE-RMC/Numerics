using System;

namespace Numerics.Data
{
    /// <summary>
    /// Enumeration of time block window options.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>

    [Serializable]
    public enum TimeBlockWindow
    {
        /// <summary>
        /// A full calendar year of 365 days from January 1st to December 31st.
        /// </summary>
        CalendarYear,

        /// <summary>
        /// A full year of 12 months from October 1st to September 30th.
        /// </summary>
        WaterYear,

        /// <summary>
        /// A custom, user specified year.
        /// </summary>
        CustomYear,

        /// <summary>
        /// Quarter of a year, a 3 month period.
        /// </summary>
        Quarter,

        /// <summary>
        /// A month of time.
        /// </summary>
        Month,       
    }
}
