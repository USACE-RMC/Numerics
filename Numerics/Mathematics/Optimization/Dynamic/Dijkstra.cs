using System;
using System.Collections.Generic;

namespace Numerics.Mathematics.Optimization
{
    /// <summary>
    /// Struct that represents an edge in a network.
    /// </summary>
    /// <remarks>
    /// Create an edge structure object. An edge contains information on the start node, end node, edge weight, and edge index.
    /// </remarks>
    /// <param name="fromNodeIndex">Node index at start of edge.</param>
    /// <param name="toNodeIndex">Node index at end of edge.</param>
    /// <param name="edgeWeight">Weight (or Cost) of the edge.</param>
    /// <param name="edgeIndex">Nonnegative index of the edge.</param>
    public struct Edge(int fromNodeIndex, int toNodeIndex, float edgeWeight, int edgeIndex)
    {
        /// <summary>
        /// Node index at start of edge.
        /// </summary>
        public int FromIndex = fromNodeIndex;
        /// <summary>
        /// Node index at end of edge.
        /// </summary>
        public int ToIndex = toNodeIndex;
        /// <summary>
        /// Weight (or Cost) of transversing the edge.
        /// </summary>
        public float Weight = edgeWeight;
        /// <summary>
        /// Nonnegative index of the edge, often used as an index to the edge source (e.g., road segment).
        /// </summary>
        public int Index = edgeIndex;
    }

    /// <summary>
    /// Dijkstra dynamic programming implementation for shortest path optimization.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// The solvers run Dijkstra's algorithm rooted at the destination, following edges backward,
    /// and return a routing table for every node: column 0 is the next node toward the
    /// destination, column 1 is the index of the edge to take, and column 2 is the cumulative
    /// cost to the destination. The destination row is (itself, -1, 0); an unreachable node's
    /// row is (-1, -1, positive infinity). Dijkstra's algorithm assumes non-negative edge
    /// weights; negative weights are not rejected, but routes computed from them are undefined.
    /// A weight of positive infinity makes an edge impassable, and a NaN weight never relaxes,
    /// severing its edge.
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// Dijkstra, E. W. (1959). A note on two problems in connexion with graphs. Numerische Mathematik, 1, 269-271.
    /// </description></item>
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm"/>
    /// </description></item>
    /// </list>
    /// </remarks>
    public static class Dijkstra
    {
        /// <summary>The result-table column holding the next node toward the destination.</summary>
        private const int NEXT_NODE = 0;

        /// <summary>The result-table column holding the index of the edge to take.</summary>
        private const int EDGE_INDEX = 1;

        /// <summary>The result-table column holding the cumulative cost to the destination.</summary>
        private const int COST = 2;

        /// <summary>
        /// Reports whether the result table records a finite-cost route from the specified node
        /// to the solved destination.
        /// </summary>
        /// <param name="resultTable">A result table produced by one of the solvers.</param>
        /// <param name="nodeIndex">The node to query.</param>
        /// <returns>True when the node can reach the destination; false when its cost is positive infinity.</returns>
        /// <exception cref="ArgumentNullException">Thrown when the result table is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the result table does not have three columns.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the node index is outside the table.</exception>
        public static bool PathExists(float[,] resultTable, int nodeIndex)
        {
            if (resultTable == null) throw new ArgumentNullException(nameof(resultTable));
            if (resultTable.GetLength(1) != 3) throw new ArgumentException("The result table must have three columns.", nameof(resultTable));
            if (nodeIndex < 0 || nodeIndex >= resultTable.GetLength(0)) throw new ArgumentOutOfRangeException(nameof(nodeIndex), $"The node index must be within [0, {resultTable.GetLength(0)}).");
            return !float.IsPositiveInfinity(resultTable[nodeIndex, COST]);
        }

        /// <summary>
        /// Walks a result table from the specified start node and returns the ordered list of
        /// edge indices leading to the destination.
        /// </summary>
        /// <param name="resultTable">A result table produced by one of the solvers.</param>
        /// <param name="startNodeIndex">The node to start from.</param>
        /// <returns>
        /// The ordered edge indices from the start node to the destination; an empty list when
        /// the start node is itself the destination; null when the destination is unreachable.
        /// </returns>
        /// <exception cref="ArgumentNullException">Thrown when the result table is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the result table does not have three columns, or does not converge to a destination.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the start node is outside the table.</exception>
        public static List<int>? GetPath(float[,] resultTable, int startNodeIndex)
        {
            return TryGetPath(resultTable, startNodeIndex, out List<int> pathEdgeIndices, out _) ? pathEdgeIndices : null;
        }

        /// <summary>
        /// Attempts to walk a result table from the specified start node, returning the ordered
        /// edge indices and the total path cost.
        /// </summary>
        /// <param name="resultTable">A result table produced by one of the solvers.</param>
        /// <param name="startNodeIndex">The node to start from.</param>
        /// <param name="pathEdgeIndices">Receives the ordered edge indices; empty when the start node is the destination or unreachable.</param>
        /// <param name="totalCost">Receives the total path cost; positive infinity when unreachable.</param>
        /// <returns>True when the destination is reachable from the start node; otherwise false.</returns>
        /// <exception cref="ArgumentNullException">Thrown when the result table is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the result table does not have three columns, or does not converge to a destination.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the start node is outside the table.</exception>
        public static bool TryGetPath(float[,] resultTable, int startNodeIndex, out List<int> pathEdgeIndices, out float totalCost)
        {
            if (resultTable == null) throw new ArgumentNullException(nameof(resultTable));
            if (resultTable.GetLength(1) != 3) throw new ArgumentException("The result table must have three columns.", nameof(resultTable));
            int rowCount = resultTable.GetLength(0);
            if (startNodeIndex < 0 || startNodeIndex >= rowCount) throw new ArgumentOutOfRangeException(nameof(startNodeIndex), $"The start node index must be within [0, {rowCount}).");

            pathEdgeIndices = new List<int>();
            totalCost = resultTable[startNodeIndex, COST];
            if (float.IsPositiveInfinity(totalCost)) return false;

            int node = startNodeIndex;
            int steps = 0;
            while (resultTable[node, EDGE_INDEX] >= 0f)
            {
                if (++steps > rowCount)
                    throw new ArgumentException("The result table does not converge to a destination; it may be inconsistent.", nameof(resultTable));
                pathEdgeIndices.Add((int)resultTable[node, EDGE_INDEX]);
                int next = (int)resultTable[node, NEXT_NODE];
                if (next < 0 || next >= rowCount)
                    throw new ArgumentException("The result table routes to a node outside the table; it may be inconsistent.", nameof(resultTable));
                node = next;
            }
            if (resultTable[node, COST] != 0f)
                throw new ArgumentException("The result table walk ended away from a destination; it may be inconsistent.", nameof(resultTable));
            return true;
        }

        /// <summary>
        /// Solves the shortest path from every node to its nearest destination in a single
        /// multi-source pass.
        /// </summary>
        /// <param name="edges">Edges, or segments, that make up the network.</param>
        /// <param name="destinationIndices">Indices of the destination nodes; at least one is required.</param>
        /// <param name="nodeCount">Optional number of nodes in the network. If not provided it will be calculated internally.</param>
        /// <returns>Lookup table of shortest paths from any given node to its nearest destination.</returns>
        /// <remarks>
        /// Costs match the multi-destination <see cref="Solve(IList{Edge}, int[], int, List{Edge}[])"/>
        /// overload exactly; the routed next node and edge can differ from it wherever two
        /// routes have exactly equal cost — whether to the same or to different destinations —
        /// where this method resolves the tie by deterministic heap order rather than
        /// destination array order. One pass over the network replaces one pass per
        /// destination. Duplicate destination indices are tolerated.
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when the edges or destination indices are null.</exception>
        /// <exception cref="ArgumentException">Thrown when the destination array is empty, the node count cannot be derived, or an edge references a node outside the network or has a negative index.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the node count is not positive, or a destination index is outside the network.</exception>
        public static float[,] SolveNearest(IList<Edge> edges, int[] destinationIndices, int nodeCount = -1)
        {
            if (edges == null) throw new ArgumentNullException(nameof(edges));
            if (destinationIndices == null) throw new ArgumentNullException(nameof(destinationIndices));
            if (destinationIndices.Length == 0) throw new ArgumentException("At least one destination index is required.", nameof(destinationIndices));

            int nNodes = ResolveNodeCount(edges, nodeCount);
            for (int i = 0; i < destinationIndices.Length; i++)
            {
                if (destinationIndices[i] < 0 || destinationIndices[i] >= nNodes)
                    throw new ArgumentOutOfRangeException(nameof(destinationIndices), $"The destination index {destinationIndices[i]} must be within [0, {nNodes}).");
            }

            CompactAdjacency adjacency = CompactAdjacency.FromEdges(edges, nNodes, groupByEndNode: true, nameof(edges));

            var next = new int[nNodes];
            var edgeIndexes = new int[nNodes];
            var dist = new float[nNodes];
            SolveNearestCore(adjacency, null, destinationIndices, next, edgeIndexes, dist, new int[nNodes], new IndexedMinHeap(nNodes));
            return WriteTable(next, edgeIndexes, dist);
        }

        /// <summary>
        /// Runs the single-pass multi-source solve: every destination seeds at cost zero (in
        /// array order, duplicates skipped) and one heap drain routes every node to its nearest
        /// destination. Buffers are fully re-initialized, so they can be reused across calls.
        /// </summary>
        /// <param name="adjacency">The incoming-edge compact adjacency.</param>
        /// <param name="weightOverride">Optional positional weight overlay; null uses the adjacency weights.</param>
        /// <param name="destinationIndices">The destination nodes.</param>
        /// <param name="next">Receives the next node toward the nearest destination per node; -1 when unreachable.</param>
        /// <param name="edgeIndexes">Receives the edge index to take per node; -1 when unreachable.</param>
        /// <param name="dist">Receives the cumulative cost per node; positive infinity when unreachable.</param>
        /// <param name="state">Scratch node states.</param>
        /// <param name="heap">The scratch heap, sized at the node count.</param>
        internal static void SolveNearestCore(in CompactAdjacency adjacency, float[]? weightOverride, int[] destinationIndices,
            int[] next, int[] edgeIndexes, float[] dist, int[] state, IndexedMinHeap heap)
        {
            int nNodes = adjacency.NodeCount;
            for (int i = 0; i < nNodes; i++)
            {
                next[i] = -1;
                edgeIndexes[i] = -1;
                dist[i] = float.PositiveInfinity;
                state[i] = 0;
            }
            heap.Clear();

            for (int d = 0; d < destinationIndices.Length; d++)
            {
                int destination = destinationIndices[d];
                if (state[destination] == 2) continue; // duplicate destination
                next[destination] = destination;
                edgeIndexes[destination] = -1;
                dist[destination] = 0f;
                heap.Add(destination, 0f);
                state[destination] = 2;
            }

            RunToExhaustion(adjacency, weightOverride, next, edgeIndexes, dist, state, heap);
        }

        /// <summary>
        /// Solves the shortest path from every node in the network of edges to the nearest of
        /// the given destinations.
        /// </summary>
        /// <param name="edges">Edges, or segments, that make up the network.</param>
        /// <param name="destinationIndices">Indices of the destination nodes.</param>
        /// <param name="nodeCount">Optional number of nodes in the network. If not provided it will be calculated internally.</param>
        /// <param name="edgesFromNodes">Optional list of incoming edges for each node in the network. If not provided or mismatched with the node count it will be calculated internally.</param>
        /// <returns>Lookup table of shortest paths from any given node.</returns>
        /// <remarks>
        /// Each destination is solved independently and the tables merge per node by strictly
        /// smaller cost in destination order, so on an exact cost tie the earlier destination in
        /// the array wins. An empty destination array returns an all-unreachable table.
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when the edges or destination indices are null.</exception>
        /// <exception cref="ArgumentException">Thrown when the node count cannot be derived, or an edge references a node outside the network or has a negative index.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the node count is not positive, or a destination index is outside the network.</exception>
        public static float[,] Solve(IList<Edge> edges, int[] destinationIndices, int nodeCount = -1, List<Edge>[]? edgesFromNodes = null)
        {
            if (edges == null) throw new ArgumentNullException(nameof(edges));
            if (destinationIndices == null) throw new ArgumentNullException(nameof(destinationIndices));

            int nNodes = ResolveNodeCount(edges, nodeCount);
            for (int i = 0; i < destinationIndices.Length; i++)
            {
                if (destinationIndices[i] < 0 || destinationIndices[i] >= nNodes)
                    throw new ArgumentOutOfRangeException(nameof(destinationIndices), $"The destination index {destinationIndices[i]} must be within [0, {nNodes}).");
            }

            CompactAdjacency adjacency = BuildIncomingAdjacency(edges, nNodes, edgesFromNodes);
            var bestNext = new int[nNodes];
            var bestEdge = new int[nNodes];
            var bestDist = new float[nNodes];
            SolveMergedCore(adjacency, null, destinationIndices,
                new int[nNodes], new int[nNodes], new float[nNodes], new int[nNodes], new IndexedMinHeap(nNodes),
                bestNext, bestEdge, bestDist);
            return WriteTable(bestNext, bestEdge, bestDist);
        }

        /// <summary>
        /// Runs one destination-rooted solve per destination in array order and merges the
        /// results per node by strictly smaller cost — the multi-destination semantics: on an
        /// exact cost tie the earlier destination wins. The accumulators start all-unreachable
        /// and the scratch buffers are reused across destinations.
        /// </summary>
        /// <param name="adjacency">The incoming-edge compact adjacency.</param>
        /// <param name="weightOverride">Optional positional weight overlay; null uses the adjacency weights.</param>
        /// <param name="destinationIndices">The destination nodes, in merge order.</param>
        /// <param name="next">Scratch per-node next-node buffer.</param>
        /// <param name="edgeIndexes">Scratch per-node edge-index buffer.</param>
        /// <param name="dist">Scratch per-node cost buffer.</param>
        /// <param name="state">Scratch per-node state buffer.</param>
        /// <param name="heap">The scratch heap, sized at the node count.</param>
        /// <param name="bestNext">Receives the merged next node per node.</param>
        /// <param name="bestEdge">Receives the merged edge index per node.</param>
        /// <param name="bestDist">Receives the merged cost per node.</param>
        internal static void SolveMergedCore(in CompactAdjacency adjacency, float[]? weightOverride, int[] destinationIndices,
            int[] next, int[] edgeIndexes, float[] dist, int[] state, IndexedMinHeap heap,
            int[] bestNext, int[] bestEdge, float[] bestDist)
        {
            int nNodes = adjacency.NodeCount;
            for (int i = 0; i < nNodes; i++)
            {
                bestNext[i] = -1;
                bestEdge[i] = -1;
                bestDist[i] = float.PositiveInfinity;
            }

            for (int d = 0; d < destinationIndices.Length; d++)
            {
                SolveCore(adjacency, weightOverride, destinationIndices[d], next, edgeIndexes, dist, state, heap);
                for (int j = 0; j < nNodes; j++)
                {
                    if (dist[j] < bestDist[j])
                    {
                        bestNext[j] = next[j];
                        bestEdge[j] = edgeIndexes[j];
                        bestDist[j] = dist[j];
                    }
                }
            }
        }

        /// <summary>
        /// Solves the shortest path from every node in the network of edges to a given destination.
        /// </summary>
        /// <param name="edges">Edges, or segments, that make up the network.</param>
        /// <param name="destinationIndex">Index of the destination node.</param>
        /// <param name="nodeCount">Optional number of nodes in the network. If not provided it will be calculated internally.</param>
        /// <param name="edgesToNodes">Optional list of incoming edges for each node in the network. If not provided or mismatched with the node count it will be calculated internally.</param>
        /// <returns>Lookup table of shortest paths from any given node.</returns>
        /// <exception cref="ArgumentNullException">Thrown when the edges are null.</exception>
        /// <exception cref="ArgumentException">Thrown when the node count cannot be derived, or an edge references a node outside the network or has a negative index.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the node count is not positive, or the destination index is outside the network.</exception>
        public static float[,] Solve(IList<Edge> edges, int destinationIndex, int nodeCount = -1, List<Edge>[]? edgesToNodes = null)
        {
            if (edges == null) throw new ArgumentNullException(nameof(edges));

            int nNodes = ResolveNodeCount(edges, nodeCount);
            if (destinationIndex < 0 || destinationIndex >= nNodes)
                throw new ArgumentOutOfRangeException(nameof(destinationIndex), $"The destination index must be within [0, {nNodes}).");

            CompactAdjacency adjacency = BuildIncomingAdjacency(edges, nNodes, edgesToNodes);

            var next = new int[nNodes];
            var edgeIndexes = new int[nNodes];
            var dist = new float[nNodes];
            var state = new int[nNodes];
            var heap = new IndexedMinHeap(nNodes);

            SolveCore(adjacency, null, destinationIndex, next, edgeIndexes, dist, state, heap);
            return WriteTable(next, edgeIndexes, dist);
        }

        /// <summary>
        /// Runs the destination-rooted solve over the compact adjacency, filling the caller's
        /// flat buffers: for each node, the next node toward the destination, the edge index to
        /// take, and the cumulative cost. Buffers are fully re-initialized, so they can be
        /// reused across calls.
        /// </summary>
        /// <param name="adjacency">The incoming-edge compact adjacency.</param>
        /// <param name="weightOverride">Optional positional weight overlay (indexed by original edge position); null uses the adjacency weights. Requires an adjacency built from an edge sequence.</param>
        /// <param name="destinationIndex">The destination node.</param>
        /// <param name="next">Receives the next node toward the destination per node; -1 when unreachable.</param>
        /// <param name="edgeIndexes">Receives the edge index to take per node; -1 when unreachable.</param>
        /// <param name="dist">Receives the cumulative cost per node; positive infinity when unreachable.</param>
        /// <param name="state">Scratch node states (0 unscanned, 1 solved, 2 in the heap).</param>
        /// <param name="heap">The scratch heap, sized at the node count.</param>
        /// <remarks>
        /// The heap holds exactly the nodes with state 2 and never a duplicate: nodes enter only
        /// through Add when their state is not 2, in-heap improvements route through
        /// DecreaseKey, and RemoveMin retires the node to state 1. The heap count therefore
        /// never exceeds the node count, independent of weight signs, which is what makes the
        /// node-count capacity exact.
        /// </remarks>
        internal static void SolveCore(in CompactAdjacency adjacency, float[]? weightOverride, int destinationIndex,
            int[] next, int[] edgeIndexes, float[] dist, int[] state, IndexedMinHeap heap)
        {
            int nNodes = adjacency.NodeCount;
            for (int i = 0; i < nNodes; i++)
            {
                next[i] = -1;
                edgeIndexes[i] = -1;
                dist[i] = float.PositiveInfinity;
                state[i] = 0;
            }
            heap.Clear();

            next[destinationIndex] = destinationIndex;
            edgeIndexes[destinationIndex] = -1;
            dist[destinationIndex] = 0f;
            heap.Add(destinationIndex, 0f);
            state[destinationIndex] = 2;

            RunToExhaustion(adjacency, weightOverride, next, edgeIndexes, dist, state, heap);
        }

        /// <summary>
        /// Drains the heap, relaxing each settled node's incoming edges — the shared main loop
        /// of the destination-rooted solvers. Seeding is the caller's responsibility.
        /// </summary>
        /// <param name="adjacency">The incoming-edge compact adjacency.</param>
        /// <param name="weightOverride">Optional positional weight overlay; null uses the adjacency weights.</param>
        /// <param name="next">The per-node next-node buffer.</param>
        /// <param name="edgeIndexes">The per-node edge-index buffer.</param>
        /// <param name="dist">The per-node cumulative-cost buffer.</param>
        /// <param name="state">The per-node state buffer.</param>
        /// <param name="heap">The seeded heap.</param>
        internal static void RunToExhaustion(in CompactAdjacency adjacency, float[]? weightOverride,
            int[] next, int[] edgeIndexes, float[] dist, int[] state, IndexedMinHeap heap)
        {
            int[] rowStart = adjacency.RowStart;
            int[] fromNodes = adjacency.FromNode;
            int[] toNodes = adjacency.ToNode;
            float[] weights = adjacency.Weight;
            int[] slotEdgeIndexes = adjacency.EdgeIndex;
            int[]? sourcePositions = adjacency.SourcePosition;

            while (heap.Count > 0)
            {
                heap.RemoveMin(out int current, out float cost);

                // Defensive only: the heap never holds a settled node (see SolveCore remarks).
                if (state[current] == 1) continue;
                state[current] = 1;

                int rowEnd = rowStart[current + 1];
                for (int k = rowStart[current]; k < rowEnd; k++)
                {
                    int from = fromNodes[k];
                    float weight = weightOverride == null ? weights[k] : weightOverride[sourcePositions![k]];
                    float newCost = cost + weight;

                    if (newCost < dist[from])
                    {
                        dist[from] = newCost;
                        if (state[from] != 2)
                        {
                            heap.Add(from, newCost);
                            state[from] = 2;
                        }
                        else
                        {
                            heap.DecreaseKey(from, newCost);
                        }
                        next[from] = toNodes[k];
                        edgeIndexes[from] = slotEdgeIndexes[k];
                    }
                }
            }
        }

        /// <summary>
        /// Resolves the working node count: -1 derives it from the largest node index the edges
        /// reference; any other value must be positive and is used as given.
        /// </summary>
        /// <param name="edges">The network edges.</param>
        /// <param name="nodeCount">The caller's node count, or -1 to derive.</param>
        /// <returns>The node count.</returns>
        /// <exception cref="ArgumentException">Thrown when the count must be derived from an empty edge list.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when an explicit count is not positive.</exception>
        private static int ResolveNodeCount(IList<Edge> edges, int nodeCount)
        {
            if (nodeCount == -1)
            {
                if (edges.Count == 0)
                    throw new ArgumentException("The node count cannot be derived from an empty edge list; provide nodeCount explicitly.", nameof(edges));
                int max = 0;
                for (int i = 0; i < edges.Count; i++)
                {
                    Edge edge = edges[i];
                    if (edge.FromIndex > max) max = edge.FromIndex;
                    if (edge.ToIndex > max) max = edge.ToIndex;
                }
                return max + 1;
            }
            if (nodeCount < 1)
                throw new ArgumentOutOfRangeException(nameof(nodeCount), "The node count must be positive, or -1 to derive it from the edges.");
            return nodeCount;
        }

        /// <summary>
        /// Builds the incoming-edge compact adjacency, honoring the documented precedence: a
        /// caller-supplied per-node list array whose length matches the node count is
        /// authoritative; otherwise the adjacency is calculated from the edges.
        /// </summary>
        /// <param name="edges">The network edges.</param>
        /// <param name="nodeCount">The working node count.</param>
        /// <param name="providedLists">The caller's per-node incoming-edge lists, or null.</param>
        /// <returns>The compact adjacency.</returns>
        private static CompactAdjacency BuildIncomingAdjacency(IList<Edge> edges, int nodeCount, List<Edge>[]? providedLists)
        {
            if (providedLists != null && providedLists.Length == nodeCount)
                return CompactAdjacency.FromIncomingLists(providedLists, nodeCount, nameof(edges));
            return CompactAdjacency.FromEdges(edges, nodeCount, groupByEndNode: true, nameof(edges));
        }

        /// <summary>
        /// Writes the flat solve buffers into the public three-column result table.
        /// </summary>
        /// <param name="next">The per-node next-node buffer.</param>
        /// <param name="edgeIndexes">The per-node edge-index buffer.</param>
        /// <param name="dist">The per-node cumulative-cost buffer.</param>
        /// <returns>The result table.</returns>
        internal static float[,] WriteTable(int[] next, int[] edgeIndexes, float[] dist)
        {
            var resultTable = new float[next.Length, 3];
            WriteTable(next, edgeIndexes, dist, resultTable);
            return resultTable;
        }

        /// <summary>
        /// Writes the flat solve buffers into a caller-supplied three-column result table.
        /// </summary>
        /// <param name="next">The per-node next-node buffer.</param>
        /// <param name="edgeIndexes">The per-node edge-index buffer.</param>
        /// <param name="dist">The per-node cumulative-cost buffer.</param>
        /// <param name="resultTable">The table to fill; its row count must equal the buffer length.</param>
        internal static void WriteTable(int[] next, int[] edgeIndexes, float[] dist, float[,] resultTable)
        {
            int nNodes = next.Length;
            for (int i = 0; i < nNodes; i++)
            {
                resultTable[i, NEXT_NODE] = next[i];
                resultTable[i, EDGE_INDEX] = edgeIndexes[i];
                resultTable[i, COST] = dist[i];
            }
        }
    }
}
