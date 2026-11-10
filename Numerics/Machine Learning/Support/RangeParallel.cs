using System;
using System.Collections.Concurrent;
using System.Threading.Tasks;

namespace Numerics.MachineLearning
{
    /// <summary>
    /// Runs independent per-index work over contiguous index ranges, in parallel only when the work
    /// is large enough to repay the scheduling cost.
    /// </summary>
    /// <remarks>
    /// Callers pass bodies that write only outputs owned by their index and perform any reduction
    /// afterwards, in index order, so results are identical whether a pass runs serially or in
    /// parallel and do not depend on the processor count.
    /// </remarks>
    internal static class RangeParallel
    {
        /// <summary>
        /// The approximate operation count below which a pass runs serially.
        /// </summary>
        internal const long WorkThreshold = 200_000;

        /// <summary>
        /// Runs <paramref name="body"/> for every index in [0, <paramref name="count"/>).
        /// </summary>
        /// <param name="count">The number of indices.</param>
        /// <param name="work">The approximate number of arithmetic operations in the pass.</param>
        /// <param name="body">The per-index body.</param>
        internal static void For(int count, long work, Action<int> body)
        {
            if (work < WorkThreshold || count < 2 || Environment.ProcessorCount == 1)
            {
                for (int i = 0; i < count; i++)
                    body(i);
                return;
            }
            int rangeSize = Math.Max(1, count / (4 * Environment.ProcessorCount));
            Parallel.ForEach(Partitioner.Create(0, count, rangeSize), range =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                    body(i);
            });
        }
    }
}
