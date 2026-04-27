namespace Numerics.Distributions
{

    /// <summary>
    /// Simple Distribution Interface.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public interface IDistribution
    {

        /// <summary>
        /// Returns the display name of the distribution type as a string.
        /// </summary>
        string DisplayName { get; }

        /// <summary>
        /// Returns the short display name of the distribution as a string.
        /// </summary>
        string ShortDisplayName { get; }
    }
}