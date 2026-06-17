using System;

namespace Numerics.Sampling
{
    /// <summary>
    /// Stores diagnostic counts and timings from the most recent pivotal bootstrap run.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class PivotalBootstrapDiagnostics
    {
        /// <summary>
        /// Gets or sets the number of raw bootstrap replicates requested or supplied.
        /// </summary>
        public int RequestedReplicates { get; set; }

        /// <summary>
        /// Gets or sets the number of raw bootstrap fits rejected before transformation.
        /// </summary>
        public int RejectedRawReplicates { get; set; }

        /// <summary>
        /// Gets or sets the number of raw bootstrap fits that failed after all retries.
        /// </summary>
        public int FailedRawReplicates { get; set; }

        /// <summary>
        /// Gets or sets the number of raw bootstrap fits accepted for transformation.
        /// </summary>
        public int AcceptedRawReplicates { get; set; }

        /// <summary>
        /// Gets or sets the number of invalid pivotal draws encountered after raw-fit acceptance.
        /// </summary>
        public int InvalidPivotalReplicates { get; set; }

        /// <summary>
        /// Gets or sets the number of pivotal draws retained after invalid-draw handling.
        /// </summary>
        public int RetainedPivotalReplicates { get; set; }

        /// <summary>
        /// Gets or sets the time spent generating and fitting raw bootstrap replicates.
        /// </summary>
        public TimeSpan ResamplingTime { get; set; }

        /// <summary>
        /// Gets or sets the time spent applying the pivotal transformation.
        /// </summary>
        public TimeSpan TransformationTime { get; set; }
    }
}
