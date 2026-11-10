using System;

namespace Numerics.Mathematics.Optimization
{

    /// <summary>
    /// The enumeration of local optimization methods for use in global optimizers.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public enum LocalMethod
    {
        /// <summary>
        /// The Adaptive Movement (Adam) optimization algorithm.
        /// </summary>
        ADAM,
        /// <summary>
        /// The Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm.
        /// </summary>
        BFGS,
        /// <summary>
        /// The Gradient Descent algorithm.
        /// </summary>
        GradientDescent,
        /// <summary>
        /// The Nelder-Mead downhill simplex algorithm
        /// </summary>
        NelderMead,
        /// <summary>
        /// The Powell optimization algorithm.
        /// </summary>
        Powell,
        /// <summary>
        /// The limited-memory BFGS algorithm for bound-constrained optimization (L-BFGS-B). Supported as the local
        /// search of <see cref="MultiStart"/> and <see cref="MLSL"/>.
        /// </summary>
        LBFGSB
    }
}
