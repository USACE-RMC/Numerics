using System.Collections.Generic;

namespace Numerics.Distributions
{
    /// <summary>
    /// An interface for estimation by the method of moments.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public interface IMomentEstimation
    {
        /// <summary>
        /// Returns an array of distribution parameters given the central moments of the sample.
        /// </summary>
        /// <param name="moments">The array of sample linear moments.</param>
        double[] ParametersFromMoments(IList<double> moments);

        /// <summary>
        /// Returns an array of central moments given the distribution parameters.
        /// </summary>
        /// <param name="parameters">The list of distribution parameters.</param>
        double[] MomentsFromParameters(IList<double> parameters);
    }
}
