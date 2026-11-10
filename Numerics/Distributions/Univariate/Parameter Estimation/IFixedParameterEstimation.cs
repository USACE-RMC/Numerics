using System.Collections.Generic;

namespace Numerics.Distributions
{

    /// <summary>
    /// Distribution parameter estimation with some parameters held at known values.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// The general form of a fit at known parameters, for example the generalized Pareto scale and shape at a known
    /// location (the threshold of a peaks-over-threshold series). A family implements the combinations of method and
    /// fixed parameters it supports and refuses the others with <see cref="System.NotSupportedException"/>. With no
    /// parameter fixed the estimate is that of <see cref="IEstimation.Estimate(IList{double}, ParameterEstimationMethod)"/>.
    /// This interface extends <see cref="IEstimation"/> rather than changing it, so existing implementations are
    /// unaffected.
    /// </para>
    /// </remarks>
    public interface IFixedParameterEstimation : IEstimation
    {

        /// <summary>
        /// Estimates and sets the parameters of the underlying distribution given a sample of observations, holding
        /// the given parameters at their values.
        /// </summary>
        /// <param name="sample">The array of sample data.</param>
        /// <param name="estimationMethod">The parameter estimation method.</param>
        /// <param name="fixedParameters">One entry per parameter, in the order of
        /// <see cref="IUnivariateDistribution.GetParameters"/>: a finite value holds that parameter fixed, and NaN
        /// marks a parameter to estimate.</param>
        void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod, IList<double> fixedParameters);
    }
}
