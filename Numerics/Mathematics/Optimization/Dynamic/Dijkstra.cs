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
    /// row is (-1, -1, positive infinity). Dijkstra's algorithm requires non-negative edge
    /// weights, and every solver throws <see cref="ArgumentOutOfRangeException"/> at entry when a
    /// weight is negative. A weight of positive infinity makes an edge impassable, and a NaN weight never relaxes,
    /// severing its edge. Every solver runs one relaxation loop in double precision: the single-precision
    /// edge weights are widened exactly, costs accumulate in double precision, and each cost is rounded
    /// to single precision once, in the result table. A cost beyond single-precision range is reported as unreachable. <see cref="ShortestDistances"/> runs the same loop on a
    /// double-precision graph in compressed sparse row form and returns the distances themselves.
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
        /// routes have exactly equal cost (whether to the same or to different destinations),
        /// where this method resolves the tie by deterministic heap order rather than
        /// destination array order. One pass over the network replaces one pass per
        /// destination. Duplicate destination indices are tolerated.
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when the edges or destination indices are null.</exception>
        /// <exception cref="ArgumentException">Thrown when the destination array is empty, the node count cannot be derived, or an edge references a node outside the network or has a negative index.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the node count is not positive, a destination index is outside the network, or an edge weight is negative.</exception>
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
            ValidateWeights(adjacency, nameof(edges));

            var next = new int[nNodes];
            var edgeIndexes = new int[nNodes];
            var dist = new double[nNodes];
            SolveNearestCore(adjacency, adjacency.Weight, destinationIndices, next, edgeIndexes, dist, new IndexedMinHeap(nNodes));
            return WriteTable(next, edgeIndexes, dist);
        }

        /// <summary>
        /// Computes the shortest distance from the nearest of several sources to every node of a directed graph in
        /// compressed sparse row (CSR) form, in double precision.
        /// </summary>
        /// <param name="nodeCount">The number of nodes; must be positive.</param>
        /// <param name="rowStart">The CSR row pointer of length <paramref name="nodeCount"/> + 1: the edges leaving node n are the
        /// slots [rowStart[n], rowStart[n + 1]). It starts at zero, never decreases, and ends at the edge count.</param>
        /// <param name="toNode">The end node of each edge slot (the CSR column indices), each within [0, <paramref name="nodeCount"/>).</param>
        /// <param name="weight">The non-negative weight of each edge slot (the CSR data). Positive infinity makes an edge impassable.</param>
        /// <param name="sources">The source nodes, each at distance zero. Duplicates are tolerated; an empty list leaves every node unreachable.</param>
        /// <param name="limit">Optional. The largest distance to compute: a node farther than this from every source is reported as
        /// unreachable. Must not be NaN. Default = positive infinity.</param>
        /// <param name="predecessors">Optional. When given, receives the node preceding each node on its shortest path, or -1 for a
        /// source or an unreachable node. Its length must be <paramref name="nodeCount"/>.</param>
        /// <param name="nearestSource">Optional. When given, receives the source each node's shortest path starts from (the node itself
        /// for a source), or -1 for an unreachable node. Its length must be <paramref name="nodeCount"/>.</param>
        /// <returns>The distance of each node from its nearest source; positive infinity when it is unreachable or beyond the limit.</returns>
        /// <exception cref="ArgumentNullException">Thrown when the row pointer, end nodes, weights or sources are null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the node count is not positive, a source or end node is outside
        /// the graph, a weight is negative or NaN, or the limit is NaN.</exception>
        /// <exception cref="ArgumentException">Thrown when the CSR arrays are inconsistent, or an output array has the wrong length.</exception>
        /// <remarks>
        /// <para>
        /// This is the multi-source form of Dijkstra's algorithm: every source is seeded at distance zero and one pass of a binary
        /// heap with decrease-key settles every node at its distance from the nearest source. It matches
        /// <c>scipy.sparse.csgraph.dijkstra(csgraph, directed=True, indices=sources, min_only=True, limit=limit,
        /// return_predecessors=True)</c>: an edge relaxes a node when the new distance is no greater than the limit and strictly less
        /// than the node's current distance, so the distances agree with SciPy. Where near-equal paths are summed in a different order,
        /// the distances can differ in the last bits, and on ties the predecessor and nearest source can differ, because the two heaps
        /// settle tied nodes in different orders.
        /// SciPy reports -9999 where this method reports -1. The edge-list solvers of this class use the same relaxation loop.
        /// </para>
        /// <para>
        /// The cost is O((V + E) log V) time for V nodes and E edges, and O(V) working memory besides the outputs.
        /// On graphs larger than the processor cache the time is dominated by memory latency, which the loop reduces by
        /// keeping each node's distance next to its heap slot and by reading the next node's edges while the current
        /// node's edges relax.
        /// </para>
        /// </remarks>
        public static double[] ShortestDistances(int nodeCount, int[] rowStart, int[] toNode, double[] weight, IList<int> sources,
            double limit = double.PositiveInfinity, int[]? predecessors = null, int[]? nearestSource = null)
        {
            if (rowStart == null) throw new ArgumentNullException(nameof(rowStart));
            if (toNode == null) throw new ArgumentNullException(nameof(toNode));
            if (weight == null) throw new ArgumentNullException(nameof(weight));
            if (sources == null) throw new ArgumentNullException(nameof(sources));
            if (nodeCount < 1) throw new ArgumentOutOfRangeException(nameof(nodeCount), "The node count must be positive.");
            if (double.IsNaN(limit)) throw new ArgumentOutOfRangeException(nameof(limit), "The limit must not be NaN.");
            if (rowStart.Length != nodeCount + 1)
                throw new ArgumentException($"The row pointer must have nodeCount + 1 = {nodeCount + 1} entries.", nameof(rowStart));
            if (toNode.Length != weight.Length)
                throw new ArgumentException($"The end node count ({toNode.Length}) must equal the weight count ({weight.Length}).", nameof(weight));
            if (rowStart[0] != 0 || rowStart[nodeCount] != toNode.Length)
                throw new ArgumentException($"The row pointer must start at 0 and end at the edge count ({toNode.Length}).", nameof(rowStart));
            for (int n = 0; n < nodeCount; n++)
            {
                if (rowStart[n + 1] < rowStart[n])
                    throw new ArgumentException($"The row pointer decreases at node {n}.", nameof(rowStart));
            }
            for (int k = 0; k < toNode.Length; k++)
            {
                if ((uint)toNode[k] >= (uint)nodeCount)
                    throw new ArgumentOutOfRangeException(nameof(toNode), $"The end node {toNode[k]} of slot {k} must be within [0, {nodeCount}).");
                ValidateWeight(weight[k], allowNaN: false, "slot", k, nameof(weight));
            }
            for (int i = 0; i < sources.Count; i++)
            {
                if ((uint)sources[i] >= (uint)nodeCount)
                    throw new ArgumentOutOfRangeException(nameof(sources), $"The source index {sources[i]} must be within [0, {nodeCount}).");
            }
            if (predecessors != null && predecessors.Length != nodeCount)
                throw new ArgumentException($"The predecessor array must have {nodeCount} entries.", nameof(predecessors));
            if (nearestSource != null && nearestSource.Length != nodeCount)
                throw new ArgumentException($"The nearest-source array must have {nodeCount} entries.", nameof(nearestSource));

            if (predecessors != null)
            {
                for (int i = 0; i < nodeCount; i++) predecessors[i] = -1;
            }
            if (nearestSource != null)
            {
                for (int i = 0; i < nodeCount; i++) nearestSource[i] = -1;
            }

            var heap = new IndexedMinHeap(nodeCount);
            for (int i = 0; i < sources.Count; i++)
            {
                int source = sources[i];
                if (heap.WeightOf(source) == 0d) continue; // duplicate source
                heap.Add(source, 0d);
                if (nearestSource != null) nearestSource[source] = source;
            }

            RunToExhaustion(rowStart, toNode, null, null, weight, limit, predecessors, null, nearestSource, heap);
            var dist = new double[nodeCount];
            CopyWeights(heap, dist);
            return dist;
        }

        /// <summary>
        /// Runs the single-pass multi-source solve: every destination seeds at cost zero (in
        /// array order, duplicates skipped) and one heap drain routes every node to its nearest
        /// destination. Buffers are fully re-initialized, so they can be reused across calls.
        /// </summary>
        /// <param name="adjacency">The incoming-edge compact adjacency.</param>
        /// <param name="slotWeights">The per-slot edge weights: the adjacency weights, or a positional overlay written in slot order.</param>
        /// <param name="destinationIndices">The destination nodes.</param>
        /// <param name="next">Receives the next node toward the nearest destination per node; -1 when unreachable.</param>
        /// <param name="edgeIndexes">Receives the edge index to take per node; -1 when unreachable.</param>
        /// <param name="dist">Receives the cumulative cost per node; positive infinity when unreachable.</param>
        /// <param name="heap">The scratch heap, sized at the node count.</param>
        internal static void SolveNearestCore(in CompactAdjacency adjacency, double[] slotWeights, int[] destinationIndices,
            int[] next, int[] edgeIndexes, double[] dist, IndexedMinHeap heap)
        {
            int nNodes = adjacency.NodeCount;
            for (int i = 0; i < nNodes; i++)
            {
                next[i] = -1;
                edgeIndexes[i] = -1;
            }
            heap.Clear();

            for (int d = 0; d < destinationIndices.Length; d++)
            {
                int destination = destinationIndices[d];
                if (heap.WeightOf(destination) == 0d) continue; // duplicate destination
                next[destination] = destination;
                edgeIndexes[destination] = -1;
                heap.Add(destination, 0d);
            }

            RunToExhaustion(adjacency, slotWeights, next, edgeIndexes, heap);
            CopyWeights(heap, dist);
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
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the node count is not positive, a destination index is outside the network, or an edge weight is negative.</exception>
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
            ValidateWeights(adjacency, nameof(edges));
            var bestNext = new int[nNodes];
            var bestEdge = new int[nNodes];
            var bestDist = new double[nNodes];
            SolveMergedCore(adjacency, adjacency.Weight, destinationIndices,
                new int[nNodes], new int[nNodes], new double[nNodes], new IndexedMinHeap(nNodes),
                bestNext, bestEdge, bestDist);
            return WriteTable(bestNext, bestEdge, bestDist);
        }

        /// <summary>
        /// Runs one destination-rooted solve per destination in array order and merges the
        /// results per node by strictly smaller cost, the multi-destination semantics: on an
        /// exact cost tie the earlier destination wins. The accumulators start all-unreachable
        /// and the scratch buffers are reused across destinations.
        /// </summary>
        /// <param name="adjacency">The incoming-edge compact adjacency.</param>
        /// <param name="slotWeights">The per-slot edge weights: the adjacency weights, or a positional overlay written in slot order.</param>
        /// <param name="destinationIndices">The destination nodes, in merge order.</param>
        /// <param name="next">Scratch per-node next-node buffer.</param>
        /// <param name="edgeIndexes">Scratch per-node edge-index buffer.</param>
        /// <param name="dist">Scratch per-node cost buffer.</param>
        /// <param name="heap">The scratch heap, sized at the node count.</param>
        /// <param name="bestNext">Receives the merged next node per node.</param>
        /// <param name="bestEdge">Receives the merged edge index per node.</param>
        /// <param name="bestDist">Receives the merged cost per node.</param>
        internal static void SolveMergedCore(in CompactAdjacency adjacency, double[] slotWeights, int[] destinationIndices,
            int[] next, int[] edgeIndexes, double[] dist, IndexedMinHeap heap,
            int[] bestNext, int[] bestEdge, double[] bestDist)
        {
            int nNodes = adjacency.NodeCount;
            for (int i = 0; i < nNodes; i++)
            {
                bestNext[i] = -1;
                bestEdge[i] = -1;
                bestDist[i] = double.PositiveInfinity;
            }

            for (int d = 0; d < destinationIndices.Length; d++)
            {
                SolveCore(adjacency, slotWeights, destinationIndices[d], next, edgeIndexes, dist, heap);
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
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the node count is not positive, the destination index is outside the network, or an edge weight is negative.</exception>
        public static float[,] Solve(IList<Edge> edges, int destinationIndex, int nodeCount = -1, List<Edge>[]? edgesToNodes = null)
        {
            if (edges == null) throw new ArgumentNullException(nameof(edges));

            int nNodes = ResolveNodeCount(edges, nodeCount);
            if (destinationIndex < 0 || destinationIndex >= nNodes)
                throw new ArgumentOutOfRangeException(nameof(destinationIndex), $"The destination index must be within [0, {nNodes}).");

            CompactAdjacency adjacency = BuildIncomingAdjacency(edges, nNodes, edgesToNodes);
            ValidateWeights(adjacency, nameof(edges));

            var next = new int[nNodes];
            var edgeIndexes = new int[nNodes];
            var dist = new double[nNodes];
            var heap = new IndexedMinHeap(nNodes);

            SolveCore(adjacency, adjacency.Weight, destinationIndex, next, edgeIndexes, dist, heap);
            return WriteTable(next, edgeIndexes, dist);
        }

        /// <summary>
        /// Runs the destination-rooted solve over the compact adjacency, filling the caller's
        /// flat buffers: for each node, the next node toward the destination, the edge index to
        /// take, and the cumulative cost. Buffers are fully re-initialized, so they can be
        /// reused across calls.
        /// </summary>
        /// <param name="adjacency">The incoming-edge compact adjacency.</param>
        /// <param name="slotWeights">The per-slot edge weights: the adjacency weights, or a positional overlay written in slot order.</param>
        /// <param name="destinationIndex">The destination node.</param>
        /// <param name="next">Receives the next node toward the destination per node; -1 when unreachable.</param>
        /// <param name="edgeIndexes">Receives the edge index to take per node; -1 when unreachable.</param>
        /// <param name="dist">Receives the cumulative cost per node; positive infinity when unreachable.</param>
        /// <param name="heap">The scratch heap, sized at the node count.</param>
        /// <remarks>
        /// The heap never holds a duplicate: a node absent from the heap enters through Add, and
        /// an improvement to a node in the heap lowers its weight in place. The heap count
        /// therefore never exceeds the node count, independent of weight signs, which is what
        /// makes the node-count capacity exact.
        /// </remarks>
        internal static void SolveCore(in CompactAdjacency adjacency, double[] slotWeights, int destinationIndex,
            int[] next, int[] edgeIndexes, double[] dist, IndexedMinHeap heap)
        {
            int nNodes = adjacency.NodeCount;
            for (int i = 0; i < nNodes; i++)
            {
                next[i] = -1;
                edgeIndexes[i] = -1;
            }
            heap.Clear();

            next[destinationIndex] = destinationIndex;
            edgeIndexes[destinationIndex] = -1;
            heap.Add(destinationIndex, 0d);

            RunToExhaustion(adjacency, slotWeights, next, edgeIndexes, heap);
            CopyWeights(heap, dist);
        }

        /// <summary>
        /// Drains the heap, relaxing each settled node's incoming edges: the shared main loop
        /// of the destination-rooted solvers. Seeding is the caller's responsibility.
        /// </summary>
        /// <param name="adjacency">The incoming-edge compact adjacency.</param>
        /// <param name="slotWeights">The per-slot edge weights: the adjacency weights, or a positional overlay written in slot order.</param>
        /// <param name="next">The per-node next-node buffer.</param>
        /// <param name="edgeIndexes">The per-node edge-index buffer.</param>
        /// <param name="heap">The seeded heap, which also holds the per-node cumulative costs.</param>
        internal static void RunToExhaustion(in CompactAdjacency adjacency, double[] slotWeights,
            int[] next, int[] edgeIndexes, IndexedMinHeap heap)
        {
            RunToExhaustion(adjacency.RowStart, adjacency.FromNode, adjacency.ToNode, adjacency.EdgeIndex, slotWeights,
                double.PositiveInfinity, next, edgeIndexes, null, heap);
        }

        /// <summary>
        /// Drains the heap over a compressed sparse row graph in double precision, the one relaxation
        /// loop behind every solver in this class. Seeding is the caller's responsibility.
        /// </summary>
        /// <param name="rowStart">The slot range of each node: node n owns slots [rowStart[n], rowStart[n + 1]).</param>
        /// <param name="relaxNode">The node each slot relaxes when its row's node settles.</param>
        /// <param name="viaNode">The node recorded in <paramref name="next"/> for each slot; null records the settled node.</param>
        /// <param name="slotEdgeIndexes">The edge index recorded in <paramref name="edgeIndexes"/> for each slot; null records nothing.</param>
        /// <param name="slotWeights">The weight of each slot.</param>
        /// <param name="limit">The largest cost relaxed; a node whose cost would exceed it stays unreached.</param>
        /// <param name="next">The per-node next-node (or predecessor) buffer, or null.</param>
        /// <param name="edgeIndexes">The per-node edge-index buffer, or null.</param>
        /// <param name="nearestSource">The per-node seed buffer, or null; a relaxed node inherits the seed of the node that reached it.</param>
        /// <param name="heap">The seeded heap. It is also the cost store: a node's cost is its last weight in the heap,
        /// positive infinity while the node is unreached.</param>
        /// <remarks>
        /// <para>
        /// A slot relaxes when its cost, the settled cost plus the slot weight, is no greater than
        /// <paramref name="limit"/> and strictly less than the node's current cost, the rule of
        /// scipy.sparse.csgraph.dijkstra. A NaN weight never relaxes and a positive infinite weight
        /// never improves a cost, so either severs its slot. The callers reject negative weights, so
        /// a removed node's cost is final and each node is removed once.
        /// </para>
        /// <para>
        /// After each removal the loop reads the first edge slot of the node now at the top of the
        /// heap, which is usually the next node settled, so that row's cache misses overlap the
        /// relaxation of the current row. The read changes no result; its value goes to an unused
        /// field so the compiler keeps the load.
        /// </para>
        /// </remarks>
        internal static void RunToExhaustion(int[] rowStart, int[] relaxNode, int[]? viaNode, int[]? slotEdgeIndexes,
            double[] slotWeights, double limit, int[]? next, int[]? edgeIndexes, int[]? nearestSource,
            IndexedMinHeap heap)
        {
            int lookAhead = 0;
            while (heap.Count > 0)
            {
                heap.RemoveMin(out int current, out double cost);
                if (heap.Count > 0)
                {
                    int ahead = rowStart[heap.MinNode];
                    if (ahead < relaxNode.Length && slotWeights[ahead] >= 0d) lookAhead += relaxNode[ahead];
                }

                int rowEnd = rowStart[current + 1];
                for (int k = rowStart[current]; k < rowEnd; k++)
                {
                    int target = relaxNode[k];
                    double newCost = cost + slotWeights[k];

                    if (newCost <= limit && newCost < heap.WeightOf(target))
                    {
                        heap.DecreaseKey(target, newCost);
                        if (next != null) next[target] = viaNode == null ? current : viaNode[k];
                        if (edgeIndexes != null) edgeIndexes[target] = slotEdgeIndexes![k];
                        if (nearestSource != null) nearestSource[target] = nearestSource[current];
                    }
                }
            }
            s_lookAheadSink = lookAhead;
        }

        /// <summary>
        /// Receives the values the relaxation loop reads ahead, so the reads are not removed as dead code. Never read.
        /// </summary>
        private static int s_lookAheadSink;

        /// <summary>
        /// Copies each node's cost, its last weight in the heap, into a per-node buffer.
        /// </summary>
        /// <param name="heap">The drained heap.</param>
        /// <param name="dist">Receives the cost of each node; its length is the node count.</param>
        private static void CopyWeights(IndexedMinHeap heap, double[] dist)
        {
            for (int i = 0; i < dist.Length; i++) dist[i] = heap.WeightOf(i);
        }

        /// <summary>
        /// Throws when a weight is negative or, unless NaN is allowed, NaN; the shared weight check of every solver.
        /// </summary>
        /// <param name="weight">The weight.</param>
        /// <param name="allowNaN">True to accept NaN, which severs its edge in the edge-list solvers.</param>
        /// <param name="location">The kind of position reported in the message, such as "slot" or "the edge at position".</param>
        /// <param name="position">The position reported in the message.</param>
        /// <param name="parameterName">The caller's parameter name for exception reporting.</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the weight is negative, or NaN and <paramref name="allowNaN"/> is false.</exception>
        internal static void ValidateWeight(double weight, bool allowNaN, string location, int position, string parameterName)
        {
            if (weight < 0d || (!allowNaN && double.IsNaN(weight)))
            {
                throw new ArgumentOutOfRangeException(parameterName, allowNaN
                    ? $"The weight of {location} {position} must be non-negative."
                    : $"The weight of {location} {position} must be non-negative and not NaN.");
            }
        }

        /// <summary>
        /// Throws when any weight of a compiled adjacency is negative. NaN is accepted.
        /// </summary>
        /// <param name="adjacency">The compiled adjacency.</param>
        /// <param name="parameterName">The caller's parameter name for exception reporting.</param>
        /// <remarks>
        /// The message names the edge's position in the original sequence when the adjacency was built from one, and
        /// otherwise the node whose caller-supplied incoming-edge list holds the edge.
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a weight is negative.</exception>
        internal static void ValidateWeights(in CompactAdjacency adjacency, string parameterName)
        {
            double[] weights = adjacency.Weight;
            int[]? sourcePosition = adjacency.SourcePosition;
            int[] rowStart = adjacency.RowStart;
            for (int n = 0; n < adjacency.NodeCount; n++)
            {
                for (int k = rowStart[n]; k < rowStart[n + 1]; k++)
                {
                    if (sourcePosition != null)
                        ValidateWeight(weights[k], allowNaN: true, "the edge at position", sourcePosition[k], parameterName);
                    else
                        ValidateWeight(weights[k], allowNaN: true, "an incoming edge of node", n, parameterName);
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
        /// <param name="dist">The per-node cumulative-cost buffer, rounded to single precision in the table.</param>
        /// <returns>The result table.</returns>
        internal static float[,] WriteTable(int[] next, int[] edgeIndexes, double[] dist)
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
        /// <param name="dist">The per-node cumulative-cost buffer, rounded to single precision in the table.</param>
        /// <param name="resultTable">The table to fill; its row count must equal the buffer length.</param>
        /// <remarks>A cost that rounds to positive infinity in single precision is written as an unreachable row.</remarks>
        internal static void WriteTable(int[] next, int[] edgeIndexes, double[] dist, float[,] resultTable)
        {
            int nNodes = next.Length;
            for (int i = 0; i < nNodes; i++)
            {
                float cost = (float)dist[i];
                // A finite double cost beyond float range rounds to infinity; report it as unreachable, (-1, -1, +inf).
                bool unreachable = float.IsPositiveInfinity(cost);
                resultTable[i, NEXT_NODE] = unreachable ? -1 : next[i];
                resultTable[i, EDGE_INDEX] = unreachable ? -1 : edgeIndexes[i];
                resultTable[i, COST] = cost;
            }
        }
    }
}
