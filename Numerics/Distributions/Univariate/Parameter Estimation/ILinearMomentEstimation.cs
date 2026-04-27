using System.Collections.Generic;

namespace Numerics.Distributions
{

    /// <summary>
    /// An interface for estimation by the method of linear moments.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public interface ILinearMomentEstimation
    {

        /// <summary>
        /// Returns an array of distribution parameters given the linear moments of the sample.
        /// </summary>
        /// <param name="moments">The array of sample linear moments.</param>
        double[] ParametersFromLinearMoments(IList<double> moments);

        /// <summary>
        /// Returns an array of linear moments given the distribution parameters.
        /// </summary>
        /// <param name="parameters">The list of distribution parameters.</param>
        double[] LinearMomentsFromParameters(IList<double> parameters);
    }
}