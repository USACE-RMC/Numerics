using System.Collections.Generic;

namespace Numerics.Distributions
{

    /// <summary>
    /// Distribution Parameter Estimation Interface.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public interface IEstimation
    {

        /// <summary>
        /// Estimates and sets the parameters of the underlying distribution given a sample of observations.
        /// </summary>
        /// <param name="sample">The array of sample data.</param>
        /// <param name="estimationMethod">The parameter estimation method.</param>
        void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod);
    }
}