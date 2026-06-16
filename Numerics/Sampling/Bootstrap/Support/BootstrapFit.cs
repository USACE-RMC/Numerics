using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.Optimization;
using System;

namespace Numerics.Sampling
{
    /// <summary>
    /// Stores a fitted parameter vector and the covariance matrix associated with that fit.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    ///     This type is used by covariance-aware bootstrap procedures, including the
    ///     pivotal bootstrap. The parameter values and covariance are cloned on input
    ///     so callers can safely reuse or mutate their source objects after construction.
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class BootstrapFit
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="BootstrapFit"/> class.
        /// </summary>
        /// <param name="parameters">The fitted parameter set.</param>
        /// <param name="covariance">The covariance matrix for <paramref name="parameters"/>.</param>
        /// <exception cref="ArgumentException">
        /// Thrown when <paramref name="parameters"/> has no values, when <paramref name="covariance"/>
        /// is not square, or when the covariance dimension does not match the parameter count.
        /// </exception>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="covariance"/> is null.</exception>
        public BootstrapFit(ParameterSet parameters, Matrix covariance)
        {
            if (parameters.Values == null || parameters.Values.Length == 0)
                throw new ArgumentException("The parameter set must contain at least one value.", nameof(parameters));
            if (covariance == null)
                throw new ArgumentNullException(nameof(covariance));
            if (!covariance.IsSquare)
                throw new ArgumentException("The covariance matrix must be square.", nameof(covariance));
            if (covariance.NumberOfRows != parameters.Values.Length)
                throw new ArgumentException("The covariance dimension must match the parameter count.", nameof(covariance));

            Parameters = parameters.Clone();
            Covariance = covariance.Clone();
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="BootstrapFit"/> class.
        /// </summary>
        /// <param name="parameters">The fitted parameter values.</param>
        /// <param name="covariance">The covariance matrix for <paramref name="parameters"/>.</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="parameters"/> or <paramref name="covariance"/> is null.</exception>
        public BootstrapFit(double[] parameters, Matrix covariance)
            : this(new ParameterSet(parameters != null ? (double[])parameters.Clone() : throw new ArgumentNullException(nameof(parameters)), double.NaN), covariance)
        {
        }

        /// <summary>
        /// Gets the fitted parameter set.
        /// </summary>
        public ParameterSet Parameters { get; }

        /// <summary>
        /// Gets the covariance matrix for <see cref="Parameters"/>.
        /// </summary>
        public Matrix Covariance { get; }

        /// <summary>
        /// Gets the number of fitted parameters.
        /// </summary>
        public int ParameterCount => Parameters.Values.Length;
    }
}
