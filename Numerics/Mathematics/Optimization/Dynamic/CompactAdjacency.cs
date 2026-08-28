using System;
using System.Collections.Generic;

namespace Numerics.Mathematics.Optimization
{
    /// <summary>
    /// A compressed sparse row (CSR) view of a network's edges, grouped by one endpoint, used by
    /// the shortest-path solvers in place of per-node edge lists: one contiguous slot range per
    /// node, with the slot's endpoints, weight, edge index, and original edge position held in
    /// parallel arrays.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// The build is a stable two-pass counting sort, so each node's slots appear in the edges'
    /// original sequence order — exactly the order per-node <c>List&lt;Edge&gt;</c> grouping
    /// produces — which keeps relaxation order, and therefore every tie outcome, identical to
    /// the list-based traversal. <see cref="SourcePosition"/> maps a slot back to its edge's
    /// position in the original sequence so a caller can overlay per-solve weights positionally
    /// without rebuilding; it is null when the view was built from caller-supplied per-node
    /// lists, which carry no positional information.
    /// </para>
    /// </remarks>
    internal readonly struct CompactAdjacency
    {
        /// <summary>The number of nodes in the network.</summary>
        internal readonly int NodeCount;

        /// <summary>The slot range per node: node n owns slots [RowStart[n], RowStart[n + 1]).</summary>
        internal readonly int[] RowStart;

        /// <summary>The slot's edge start node (<see cref="Edge.FromIndex"/>).</summary>
        internal readonly int[] FromNode;

        /// <summary>The slot's edge end node (<see cref="Edge.ToIndex"/>).</summary>
        internal readonly int[] ToNode;

        /// <summary>The slot's edge weight.</summary>
        internal readonly float[] Weight;

        /// <summary>The slot's edge index (<see cref="Edge.Index"/>).</summary>
        internal readonly int[] EdgeIndex;

        /// <summary>
        /// The slot's edge position in the original edge sequence, for positional weight
        /// overlays; null when the view was built from caller-supplied per-node lists.
        /// </summary>
        internal readonly int[]? SourcePosition;

        /// <summary>
        /// Creates the view over the given arrays.
        /// </summary>
        /// <param name="nodeCount">The number of nodes.</param>
        /// <param name="rowStart">The per-node slot ranges.</param>
        /// <param name="fromNode">The per-slot edge start nodes.</param>
        /// <param name="toNode">The per-slot edge end nodes.</param>
        /// <param name="weight">The per-slot edge weights.</param>
        /// <param name="edgeIndex">The per-slot edge indices.</param>
        /// <param name="sourcePosition">The per-slot original edge positions, or null when unavailable.</param>
        private CompactAdjacency(int nodeCount, int[] rowStart, int[] fromNode, int[] toNode,
            float[] weight, int[] edgeIndex, int[]? sourcePosition)
        {
            NodeCount = nodeCount;
            RowStart = rowStart;
            FromNode = fromNode;
            ToNode = toNode;
            Weight = weight;
            EdgeIndex = edgeIndex;
            SourcePosition = sourcePosition;
        }

        /// <summary>
        /// Builds the view from an edge sequence, grouping by the end node (incoming edges, the
        /// orientation the destination-rooted solver walks) or by the start node (outgoing
        /// edges, the orientation forward routing walks).
        /// </summary>
        /// <param name="edges">The edges that make up the network.</param>
        /// <param name="nodeCount">The number of nodes; every referenced node index must lie in [0, nodeCount).</param>
        /// <param name="groupByEndNode">True to group slots by <see cref="Edge.ToIndex"/>; false to group by <see cref="Edge.FromIndex"/>.</param>
        /// <param name="parameterName">The caller's parameter name for exception reporting.</param>
        /// <returns>The compact view.</returns>
        /// <exception cref="ArgumentException">Thrown when an edge references a node index outside [0, nodeCount).</exception>
        internal static CompactAdjacency FromEdges(IList<Edge> edges, int nodeCount, bool groupByEndNode, string parameterName)
        {
            int edgeCount = edges.Count;
            var rowStart = new int[nodeCount + 1];

            for (int i = 0; i < edgeCount; i++)
            {
                Edge edge = edges[i];
                if (edge.FromIndex < 0 || edge.FromIndex >= nodeCount || edge.ToIndex < 0 || edge.ToIndex >= nodeCount)
                {
                    throw new ArgumentException($"The edge at position {i} references node index " +
                        $"{(edge.FromIndex < 0 || edge.FromIndex >= nodeCount ? edge.FromIndex : edge.ToIndex)}, " +
                        $"outside the node range [0, {nodeCount}).", parameterName);
                }
                rowStart[(groupByEndNode ? edge.ToIndex : edge.FromIndex) + 1]++;
            }
            for (int n = 0; n < nodeCount; n++) rowStart[n + 1] += rowStart[n];

            var cursor = new int[nodeCount];
            var fromNode = new int[edgeCount];
            var toNode = new int[edgeCount];
            var weight = new float[edgeCount];
            var edgeIndex = new int[edgeCount];
            var sourcePosition = new int[edgeCount];

            for (int i = 0; i < edgeCount; i++)
            {
                Edge edge = edges[i];
                int bucket = groupByEndNode ? edge.ToIndex : edge.FromIndex;
                int slot = rowStart[bucket] + cursor[bucket];
                cursor[bucket]++;
                fromNode[slot] = edge.FromIndex;
                toNode[slot] = edge.ToIndex;
                weight[slot] = edge.Weight;
                edgeIndex[slot] = edge.Index;
                sourcePosition[slot] = i;
            }

            return new CompactAdjacency(nodeCount, rowStart, fromNode, toNode, weight, edgeIndex, sourcePosition);
        }

        /// <summary>
        /// Builds the view from caller-supplied per-node incoming-edge lists, preserving each
        /// list's order. The lists are authoritative when supplied (the documented solver
        /// precedence), so their contents are read as-is; only out-of-range node indices carried
        /// by the listed edges throw. An edge placed in the wrong per-node bucket is not
        /// detected and silently corrupts routing — bucket consistency is the caller's
        /// responsibility.
        /// </summary>
        /// <param name="edgesToNodes">The incoming edges for each node; a null entry means the node has none.</param>
        /// <param name="nodeCount">The number of nodes; the array length must equal it.</param>
        /// <param name="parameterName">The caller's parameter name for exception reporting.</param>
        /// <returns>The compact view, without source positions.</returns>
        /// <exception cref="ArgumentException">Thrown when a listed edge references a node index outside [0, nodeCount).</exception>
        internal static CompactAdjacency FromIncomingLists(List<Edge>[] edgesToNodes, int nodeCount, string parameterName)
        {
            var rowStart = new int[nodeCount + 1];
            int edgeCount = 0;
            for (int n = 0; n < nodeCount; n++)
            {
                int count = edgesToNodes[n]?.Count ?? 0;
                rowStart[n + 1] = rowStart[n] + count;
                edgeCount += count;
            }

            var fromNode = new int[edgeCount];
            var toNode = new int[edgeCount];
            var weight = new float[edgeCount];
            var edgeIndex = new int[edgeCount];

            for (int n = 0; n < nodeCount; n++)
            {
                List<Edge>? list = edgesToNodes[n];
                if (list == null) continue;
                int slot = rowStart[n];
                for (int j = 0; j < list.Count; j++)
                {
                    Edge edge = list[j];
                    if (edge.FromIndex < 0 || edge.FromIndex >= nodeCount || edge.ToIndex < 0 || edge.ToIndex >= nodeCount)
                    {
                        throw new ArgumentException($"The incoming-edge list for node {n} contains an edge referencing node index " +
                            $"{(edge.FromIndex < 0 || edge.FromIndex >= nodeCount ? edge.FromIndex : edge.ToIndex)}, " +
                            $"outside the node range [0, {nodeCount}).", parameterName);
                    }
                    fromNode[slot] = edge.FromIndex;
                    toNode[slot] = edge.ToIndex;
                    weight[slot] = edge.Weight;
                    edgeIndex[slot] = edge.Index;
                    slot++;
                }
            }

            return new CompactAdjacency(nodeCount, rowStart, fromNode, toNode, weight, edgeIndex, null);
        }
    }
}
