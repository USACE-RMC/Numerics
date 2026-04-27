using System;

namespace Numerics.Mathematics.Optimization
{
    /// <summary>
    /// Enumeration of optimization statuses. 
    /// </summary>
    ///  <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public enum OptimizationStatus
    {
        /// <summary>
        /// Optimization has not been performed yet. 
        /// </summary>
        None, 
        
        /// <summary>
        /// The optimization method ended successfully.
        /// </summary>
        Success,

        /// <summary>
        /// The optimization method was stopped because the maximum number of iterations was reached. 
        /// </summary>
        MaximumIterationsReached,

        /// <summary>
        /// The optimization method was stopped because the maximum number of objective function evaluations was reached. 
        /// </summary>
        MaximumFunctionEvaluationsReached,

        /// <summary>
        /// The optimization method was stopped due to internal failure. 
        /// </summary>
        Failure,
    }
}
