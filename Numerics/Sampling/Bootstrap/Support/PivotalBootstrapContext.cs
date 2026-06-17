using System;

namespace Numerics.Sampling
{
    /// <summary>
    /// Provides accepted raw bootstrap fits to a pivotal bootstrap link factory.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    ///     Link factories can use this context to select fixed links, such as
    ///     <c>LogLink</c> for positive scale parameters, or to fit data-adaptive links,
    ///     such as <c>YeoJohnsonLink</c>, from the accepted raw bootstrap ensemble.
    /// </para>
    /// </remarks>
    public sealed class PivotalBootstrapContext
    {
        private readonly double[,] _rawParameterValues;

        /// <summary>
        /// Initializes a new instance of the <see cref="PivotalBootstrapContext"/> class.
        /// </summary>
        /// <param name="parentFit">The original parent fit.</param>
        /// <param name="rawBootstrapFits">The accepted raw bootstrap fits.</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="parentFit"/> or <paramref name="rawBootstrapFits"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when a raw fit parameter count differs from the parent fit.</exception>
        public PivotalBootstrapContext(BootstrapFit parentFit, BootstrapFit[] rawBootstrapFits)
        {
            ParentFit = parentFit ?? throw new ArgumentNullException(nameof(parentFit));
            RawBootstrapFits = rawBootstrapFits ?? throw new ArgumentNullException(nameof(rawBootstrapFits));

            int p = parentFit.ParameterCount;
            _rawParameterValues = new double[rawBootstrapFits.Length, p];
            for (int i = 0; i < rawBootstrapFits.Length; i++)
            {
                if (rawBootstrapFits[i].ParameterCount != p)
                    throw new ArgumentException("Every raw bootstrap fit must have the same parameter count as the parent fit.", nameof(rawBootstrapFits));

                for (int j = 0; j < p; j++)
                    _rawParameterValues[i, j] = rawBootstrapFits[i].Parameters.Values[j];
            }
        }

        /// <summary>
        /// Gets the original parent fit.
        /// </summary>
        public BootstrapFit ParentFit { get; }

        /// <summary>
        /// Gets the accepted raw bootstrap fits.
        /// </summary>
        public BootstrapFit[] RawBootstrapFits { get; }

        /// <summary>
        /// Gets the number of fitted parameters.
        /// </summary>
        public int ParameterCount => ParentFit.ParameterCount;

        /// <summary>
        /// Gets a copy of the accepted raw bootstrap values for one parameter.
        /// </summary>
        /// <param name="parameterIndex">The zero-based parameter index.</param>
        /// <returns>The accepted raw bootstrap values for <paramref name="parameterIndex"/>.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="parameterIndex"/> is outside the parameter range.</exception>
        public double[] GetRawParameterValues(int parameterIndex)
        {
            if (parameterIndex < 0 || parameterIndex >= ParameterCount)
                throw new ArgumentOutOfRangeException(nameof(parameterIndex));

            var values = new double[RawBootstrapFits.Length];
            for (int i = 0; i < values.Length; i++)
                values[i] = _rawParameterValues[i, parameterIndex];
            return values;
        }
    }
}
