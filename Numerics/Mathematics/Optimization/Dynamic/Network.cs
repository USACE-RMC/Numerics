using System;
using System.Collections.Generic;

namespace Numerics.Mathematics.Optimization
{
    /// <summary>
    /// A network of edges used for shortest path optimization applications.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// The network compiles its topology once at construction — the node count, the incoming and
    /// outgoing adjacency, and the destination set are fixed for the instance's lifetime — so
    /// repeated solves pay only the solve itself, and custom-weight solves overlay a positional
    /// weight vector with no rebuild. The solve and path methods that return their results
    /// allocate them per call and are safe for concurrent use; the overloads that write into a
    /// caller-supplied table reuse instance scratch buffers and are not thread safe. Weights follow the
    /// <see cref="Dijkstra"/> conventions: non-negative weights are the correctness
    /// precondition, positive infinity is impassable, and NaN severs its edge.
    /// </para>
    /// </remarks>
    public class Network
    {
        /// <summary>The outgoing edges per node, in original edge order.</summary>
        private readonly List<Edge>[] _outgoingEdges;

        /// <summary>The incoming edges per node, in original edge order.</summary>
        private readonly List<Edge>[] _incomingEdges;

        /// <summary>The number of nodes, one past the largest referenced node index.</summary>
        private readonly int _nodeCount;

        /// <summary>The network's destination nodes.</summary>
        private readonly int[] _destinationIndices;

        /// <summary>The network's edges, in caller order.</summary>
        private readonly Edge[] _edges;

        /// <summary>The compiled incoming adjacency (grouped by end node) the solvers walk.</summary>
        private readonly CompactAdjacency _incomingAdjacency;

        /// <summary>The compiled outgoing adjacency (grouped by start node) forward routing walks.</summary>
        private readonly CompactAdjacency _outgoingAdjacency;

        /// <summary>Marks the destination nodes for O(1) membership tests.</summary>
        private readonly bool[] _isDestination;

        /// <summary>Lazily allocated scratch for the table-reuse solves (not thread safe).</summary>
        private int[]? _scratchNext;

        /// <summary>Scratch per-node edge indices for the table-reuse solves.</summary>
        private int[]? _scratchEdge;

        /// <summary>Scratch per-node costs for the table-reuse solves.</summary>
        private float[]? _scratchDist;

        /// <summary>Scratch per-node states for the table-reuse solves.</summary>
        private int[]? _scratchState;

        /// <summary>Scratch heap for the table-reuse solves.</summary>
        private IndexedMinHeap? _scratchHeap;

        /// <summary>Scratch merged next nodes for the table-reuse multi-destination solve.</summary>
        private int[]? _scratchBestNext;

        /// <summary>Scratch merged edge indices for the table-reuse multi-destination solve.</summary>
        private int[]? _scratchBestEdge;

        /// <summary>Scratch merged costs for the table-reuse multi-destination solve.</summary>
        private float[]? _scratchBestDist;

        /// <summary>
        /// The destination node indices for shortest path computation.
        /// </summary>
        public int[] DestinationIndices { get => _destinationIndices; }

        /// <summary>
        /// The incoming edges for each node, indexed by node index. Nodes with no incoming
        /// edges hold an empty list.
        /// </summary>
        public List<Edge>[] IncomingEdges { get => _incomingEdges; }

        /// <summary>
        /// The outgoing edges for each node, indexed by node index. Nodes with no outgoing
        /// edges hold an empty list.
        /// </summary>
        public List<Edge>[] OutgoingEdges { get => _outgoingEdges; }

        /// <summary>
        /// The number of nodes in the network, one past the largest node index the edges
        /// reference. Caller-supplied result tables must have this many rows.
        /// </summary>
        public int NodeCount { get => _nodeCount; }

        /// <summary>
        /// Creates a new network from the specified edges and destination indices.
        /// </summary>
        /// <param name="edges">The edges that define the network. The array is copied.</param>
        /// <param name="destinationIndices">The destination node indices. The array is copied.</param>
        /// <exception cref="ArgumentNullException">Thrown when either array is null.</exception>
        /// <exception cref="ArgumentException">Thrown when either array is empty, or an edge references a negative node index.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a destination index is outside the network.</exception>
        public Network(Edge[] edges, int[] destinationIndices)
        {
            if (edges == null) throw new ArgumentNullException(nameof(edges));
            if (edges.Length == 0) throw new ArgumentException("At least one edge is required.", nameof(edges));
            if (destinationIndices == null) throw new ArgumentNullException(nameof(destinationIndices));
            if (destinationIndices.Length == 0) throw new ArgumentException("At least one destination index is required.", nameof(destinationIndices));

            _edges = new Edge[edges.Length];
            Array.Copy(edges, _edges, edges.Length);
            _destinationIndices = new int[destinationIndices.Length];
            Array.Copy(destinationIndices, _destinationIndices, destinationIndices.Length);

            int max = 0;
            for (int i = 0; i < _edges.Length; i++)
            {
                if (_edges[i].FromIndex > max) max = _edges[i].FromIndex;
                if (_edges[i].ToIndex > max) max = _edges[i].ToIndex;
            }
            _nodeCount = max + 1;

            // The compiled builds also range-check every referenced node index.
            _incomingAdjacency = CompactAdjacency.FromEdges(_edges, _nodeCount, groupByEndNode: true, nameof(edges));
            _outgoingAdjacency = CompactAdjacency.FromEdges(_edges, _nodeCount, groupByEndNode: false, nameof(edges));

            _isDestination = new bool[_nodeCount];
            for (int d = 0; d < _destinationIndices.Length; d++)
            {
                int destination = _destinationIndices[d];
                if (destination < 0 || destination >= _nodeCount)
                    throw new ArgumentOutOfRangeException(nameof(destinationIndices), $"The destination index {destination} must be within [0, {_nodeCount}).");
                _isDestination[destination] = true;
            }

            _incomingEdges = MaterializeLists(_incomingAdjacency);
            _outgoingEdges = MaterializeLists(_outgoingAdjacency);
        }

        /// <summary>
        /// Solves the shortest path from all nodes to the specified destination.
        /// </summary>
        /// <param name="destinationIndex">The destination node index.</param>
        /// <returns>A result table with the next node, edge index, and cumulative weight for each node.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the destination index is outside the network.</exception>
        public float[,] Solve(int destinationIndex)
        {
            if (destinationIndex < 0 || destinationIndex >= _nodeCount)
                throw new ArgumentOutOfRangeException(nameof(destinationIndex), $"The destination index must be within [0, {_nodeCount}).");

            var next = new int[_nodeCount];
            var edgeIndexes = new int[_nodeCount];
            var dist = new float[_nodeCount];
            Dijkstra.SolveCore(_incomingAdjacency, null, destinationIndex, next, edgeIndexes, dist, new int[_nodeCount], new IndexedMinHeap(_nodeCount));
            return Dijkstra.WriteTable(next, edgeIndexes, dist);
        }

        /// <summary>
        /// Solves the shortest path from all nodes to the nearest of the specified destinations.
        /// </summary>
        /// <param name="destinationIndices">The destination node indices.</param>
        /// <returns>A result table with the next node, edge index, and cumulative weight for each node.</returns>
        /// <remarks>
        /// Each destination is solved independently and the tables merge per node by strictly
        /// smaller cost in destination order, so on an exact cost tie the earlier destination
        /// wins — the same semantics as the static multi-destination solver. An empty
        /// destination array returns an all-unreachable table.
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when the destination indices are null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a destination index is outside the network.</exception>
        public float[,] Solve(int[] destinationIndices)
        {
            ValidateDestinations(destinationIndices);
            var bestNext = new int[_nodeCount];
            var bestEdge = new int[_nodeCount];
            var bestDist = new float[_nodeCount];
            Dijkstra.SolveMergedCore(_incomingAdjacency, null, destinationIndices,
                new int[_nodeCount], new int[_nodeCount], new float[_nodeCount], new int[_nodeCount], new IndexedMinHeap(_nodeCount),
                bestNext, bestEdge, bestDist);
            return Dijkstra.WriteTable(bestNext, bestEdge, bestDist);
        }

        /// <summary>
        /// Solves the shortest path to the network's destinations using custom edge weights.
        /// </summary>
        /// <param name="edgeWeights">Custom weights, one per edge, positional with the constructor's edge array.</param>
        /// <returns>A result table with the next node, edge index, and cumulative weight for each node.</returns>
        /// <remarks>
        /// The weights overlay the compiled topology positionally — nothing is rebuilt. The
        /// destinations are the constructor's, merged in destination order exactly like
        /// <see cref="Solve(int[])"/>.
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when the weights are null.</exception>
        /// <exception cref="ArgumentException">Thrown when the weight count does not equal the edge count.</exception>
        public float[,] Solve(float[] edgeWeights)
        {
            ValidateWeights(edgeWeights);
            var bestNext = new int[_nodeCount];
            var bestEdge = new int[_nodeCount];
            var bestDist = new float[_nodeCount];
            Dijkstra.SolveMergedCore(_incomingAdjacency, edgeWeights, _destinationIndices,
                new int[_nodeCount], new int[_nodeCount], new float[_nodeCount], new int[_nodeCount], new IndexedMinHeap(_nodeCount),
                bestNext, bestEdge, bestDist);
            return Dijkstra.WriteTable(bestNext, bestEdge, bestDist);
        }

        /// <summary>
        /// Solves the shortest path from every node to its nearest network destination in a
        /// single multi-source pass.
        /// </summary>
        /// <returns>A result table with the next node, edge index, and cumulative weight for each node.</returns>
        /// <remarks>
        /// Costs match <see cref="Solve(int[])"/> over the network's destinations exactly; the
        /// routed next node and edge can differ from it wherever two routes have exactly equal
        /// cost — whether to the same or to different destinations — where this method resolves
        /// the tie by deterministic heap order rather than destination array order. One pass
        /// replaces one pass per destination.
        /// </remarks>
        public float[,] SolveNearest()
        {
            var next = new int[_nodeCount];
            var edgeIndexes = new int[_nodeCount];
            var dist = new float[_nodeCount];
            Dijkstra.SolveNearestCore(_incomingAdjacency, null, _destinationIndices, next, edgeIndexes, dist, new int[_nodeCount], new IndexedMinHeap(_nodeCount));
            return Dijkstra.WriteTable(next, edgeIndexes, dist);
        }

        /// <summary>
        /// Solves the shortest path to the network's destinations using custom edge weights,
        /// writing the results into a caller-supplied table. Reuses instance scratch buffers —
        /// zero allocation per call at steady state — and is therefore not thread safe.
        /// </summary>
        /// <param name="edgeWeights">Custom weights, one per edge, positional with the constructor's edge array.</param>
        /// <param name="resultTable">The table to fill; its dimensions must be [<see cref="NodeCount"/>, 3].</param>
        /// <exception cref="ArgumentNullException">Thrown when the weights or table are null.</exception>
        /// <exception cref="ArgumentException">Thrown when the weight count or table dimensions are wrong.</exception>
        public void Solve(float[] edgeWeights, float[,] resultTable)
        {
            ValidateWeights(edgeWeights);
            ValidateResultTable(resultTable);
            EnsureScratch();
            Dijkstra.SolveMergedCore(_incomingAdjacency, edgeWeights, _destinationIndices,
                _scratchNext!, _scratchEdge!, _scratchDist!, _scratchState!, _scratchHeap!,
                _scratchBestNext!, _scratchBestEdge!, _scratchBestDist!);
            Dijkstra.WriteTable(_scratchBestNext!, _scratchBestEdge!, _scratchBestDist!, resultTable);
        }

        /// <summary>
        /// Solves the shortest path from every node to its nearest network destination using
        /// custom edge weights, writing the results into a caller-supplied table. Reuses
        /// instance scratch buffers — zero allocation per call at steady state — and is
        /// therefore not thread safe.
        /// </summary>
        /// <param name="edgeWeights">Custom weights, one per edge, positional with the constructor's edge array.</param>
        /// <param name="resultTable">The table to fill; its dimensions must be [<see cref="NodeCount"/>, 3].</param>
        /// <exception cref="ArgumentNullException">Thrown when the weights or table are null.</exception>
        /// <exception cref="ArgumentException">Thrown when the weight count or table dimensions are wrong.</exception>
        public void SolveNearest(float[] edgeWeights, float[,] resultTable)
        {
            ValidateWeights(edgeWeights);
            ValidateResultTable(resultTable);
            EnsureScratch();
            Dijkstra.SolveNearestCore(_incomingAdjacency, edgeWeights, _destinationIndices,
                _scratchNext!, _scratchEdge!, _scratchDist!, _scratchState!, _scratchHeap!);
            Dijkstra.WriteTable(_scratchNext!, _scratchEdge!, _scratchDist!, resultTable);
        }

        /// <summary>
        /// Materializes the per-node edge lists the public properties expose from a compiled
        /// adjacency, preserving slot order and giving every node a list.
        /// </summary>
        /// <param name="adjacency">The compiled adjacency.</param>
        /// <returns>The per-node edge lists.</returns>
        private static List<Edge>[] MaterializeLists(in CompactAdjacency adjacency)
        {
            var lists = new List<Edge>[adjacency.NodeCount];
            for (int n = 0; n < adjacency.NodeCount; n++)
            {
                int start = adjacency.RowStart[n];
                int end = adjacency.RowStart[n + 1];
                var list = new List<Edge>(end - start);
                for (int k = start; k < end; k++)
                {
                    list.Add(new Edge(adjacency.FromNode[k], adjacency.ToNode[k], adjacency.Weight[k], adjacency.EdgeIndex[k]));
                }
                lists[n] = list;
            }
            return lists;
        }

        /// <summary>
        /// Validates a caller-supplied destination array against the network.
        /// </summary>
        /// <param name="destinationIndices">The destination node indices.</param>
        /// <exception cref="ArgumentNullException">Thrown when the array is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when an index is outside the network.</exception>
        private void ValidateDestinations(int[] destinationIndices)
        {
            if (destinationIndices == null) throw new ArgumentNullException(nameof(destinationIndices));
            for (int i = 0; i < destinationIndices.Length; i++)
            {
                if (destinationIndices[i] < 0 || destinationIndices[i] >= _nodeCount)
                    throw new ArgumentOutOfRangeException(nameof(destinationIndices), $"The destination index {destinationIndices[i]} must be within [0, {_nodeCount}).");
            }
        }

        /// <summary>
        /// Validates a positional custom-weight vector against the network's edge count.
        /// </summary>
        /// <param name="edgeWeights">The custom weights.</param>
        /// <exception cref="ArgumentNullException">Thrown when the weights are null.</exception>
        /// <exception cref="ArgumentException">Thrown when the weight count does not equal the edge count.</exception>
        private void ValidateWeights(float[] edgeWeights)
        {
            if (edgeWeights == null) throw new ArgumentNullException(nameof(edgeWeights));
            if (edgeWeights.Length != _edges.Length)
                throw new ArgumentException($"The weight count ({edgeWeights.Length}) must equal the edge count ({_edges.Length}).", nameof(edgeWeights));
        }

        /// <summary>
        /// Validates a caller-supplied result table's dimensions.
        /// </summary>
        /// <param name="resultTable">The table to fill.</param>
        /// <exception cref="ArgumentNullException">Thrown when the table is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the dimensions are not [NodeCount, 3].</exception>
        private void ValidateResultTable(float[,] resultTable)
        {
            if (resultTable == null) throw new ArgumentNullException(nameof(resultTable));
            if (resultTable.GetLength(0) != _nodeCount || resultTable.GetLength(1) != 3)
                throw new ArgumentException($"The result table must be [{_nodeCount}, 3].", nameof(resultTable));
        }

        /// <summary>
        /// Allocates the reuse-solve scratch buffers on first use.
        /// </summary>
        private void EnsureScratch()
        {
            if (_scratchNext != null) return;
            _scratchNext = new int[_nodeCount];
            _scratchEdge = new int[_nodeCount];
            _scratchDist = new float[_nodeCount];
            _scratchState = new int[_nodeCount];
            _scratchHeap = new IndexedMinHeap(_nodeCount);
            _scratchBestNext = new int[_nodeCount];
            _scratchBestEdge = new int[_nodeCount];
            _scratchBestDist = new float[_nodeCount];
        }

        /// <summary>
        /// Finds an alternative path avoiding the specified edges.
        /// </summary>
        /// <param name="edgesToRemove">Edge indices to exclude from the path. The array is not modified; every edge bearing a listed index is excluded.</param>
        /// <param name="startNodeIndex">The starting node index.</param>
        /// <returns>
        /// The ordered edge indices from the start node to its nearest network destination
        /// avoiding the excluded edges; an empty list when the start node is itself a
        /// destination; null if no path exists.
        /// </returns>
        /// <exception cref="ArgumentNullException">Thrown when the edge indices are null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the start node is outside the network.</exception>
        public List<int>? GetPath(int[] edgesToRemove, int startNodeIndex)
        {
            if (edgesToRemove == null) throw new ArgumentNullException(nameof(edgesToRemove));
            if (startNodeIndex < 0 || startNodeIndex >= _nodeCount)
                throw new ArgumentOutOfRangeException(nameof(startNodeIndex), $"The start node index must be within [0, {_nodeCount}).");

            var removed = new HashSet<int>();
            for (int i = 0; i < edgesToRemove.Length; i++) removed.Add(edgesToRemove[i]);

            return FindDetourPath(removed, startNodeIndex);
        }

        /// <summary>
        /// Finds an alternative path avoiding the specified edges, using a pre-computed results
        /// table to skip the solve when the recorded route is unaffected.
        /// </summary>
        /// <param name="edgesToRemove">Edge indices to exclude from the path. The array is not modified; every edge bearing a listed index is excluded.</param>
        /// <param name="startNodeIndex">The starting node index.</param>
        /// <param name="existingResultsTable">A result table previously solved on this network toward its destinations, without exclusions.</param>
        /// <returns>
        /// The ordered edge indices from the start node to its nearest network destination
        /// avoiding the excluded edges; an empty list when the start node is itself a
        /// destination or when no path exists.
        /// </returns>
        /// <remarks>
        /// When the table's recorded route from the start node avoids every excluded edge it is
        /// returned directly — exclusions only remove paths, so a surviving unexcluded optimum
        /// stays optimal. Otherwise the path is re-solved with the exclusions applied.
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when the edge indices or table are null.</exception>
        /// <exception cref="ArgumentException">Thrown when the table dimensions are not [<see cref="NodeCount"/>, 3].</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the start node is outside the network.</exception>
        public List<int>? GetPath(int[] edgesToRemove, int startNodeIndex, float[,] existingResultsTable)
        {
            if (edgesToRemove == null) throw new ArgumentNullException(nameof(edgesToRemove));
            ValidateResultTable(existingResultsTable);
            if (startNodeIndex < 0 || startNodeIndex >= _nodeCount)
                throw new ArgumentOutOfRangeException(nameof(startNodeIndex), $"The start node index must be within [0, {_nodeCount}).");

            // Exclusions only shrink reachability: unreachable without them means unreachable with them.
            if (float.IsPositiveInfinity(existingResultsTable[startNodeIndex, 2])) return new List<int>();

            var removed = new HashSet<int>();
            for (int i = 0; i < edgesToRemove.Length; i++) removed.Add(edgesToRemove[i]);

            List<int>? recorded = Dijkstra.GetPath(existingResultsTable, startNodeIndex);
            if (recorded != null)
            {
                bool blocked = false;
                for (int i = 0; i < recorded.Count; i++)
                {
                    if (removed.Contains(recorded[i]))
                    {
                        blocked = true;
                        break;
                    }
                }
                if (!blocked) return recorded;
            }

            return FindDetourPath(removed, startNodeIndex) ?? new List<int>();
        }

        /// <summary>
        /// Runs a forward Dijkstra search from the start node over the outgoing adjacency,
        /// skipping excluded edges, stopping at the first settled destination (the nearest one),
        /// and reconstructing the edge-index path.
        /// </summary>
        /// <param name="removed">The excluded edge indices.</param>
        /// <param name="startNodeIndex">The starting node index.</param>
        /// <returns>The ordered edge indices to the nearest destination; an empty list when the start node is a destination; null when every destination is unreachable.</returns>
        private List<int>? FindDetourPath(HashSet<int> removed, int startNodeIndex)
        {
            if (_isDestination[startNodeIndex]) return new List<int>();

            var dist = new float[_nodeCount];
            var state = new int[_nodeCount];
            var previousSlot = new int[_nodeCount];
            for (int i = 0; i < _nodeCount; i++)
            {
                dist[i] = float.PositiveInfinity;
                state[i] = 0;
                previousSlot[i] = -1;
            }
            var heap = new IndexedMinHeap(_nodeCount);

            dist[startNodeIndex] = 0f;
            heap.Add(startNodeIndex, 0f);
            state[startNodeIndex] = 2;

            int[] rowStart = _outgoingAdjacency.RowStart;
            int[] toNodes = _outgoingAdjacency.ToNode;
            float[] weights = _outgoingAdjacency.Weight;
            int[] edgeIndexes = _outgoingAdjacency.EdgeIndex;

            int reachedDestination = -1;
            while (heap.Count > 0)
            {
                heap.RemoveMin(out int current, out float cost);
                if (state[current] == 1) continue;
                state[current] = 1;

                if (_isDestination[current])
                {
                    reachedDestination = current;
                    break;
                }

                int rowEnd = rowStart[current + 1];
                for (int k = rowStart[current]; k < rowEnd; k++)
                {
                    if (removed.Contains(edgeIndexes[k])) continue;
                    int to = toNodes[k];
                    float newCost = cost + weights[k];
                    if (newCost < dist[to])
                    {
                        dist[to] = newCost;
                        if (state[to] != 2)
                        {
                            heap.Add(to, newCost);
                            state[to] = 2;
                        }
                        else
                        {
                            heap.DecreaseKey(to, newCost);
                        }
                        previousSlot[to] = k;
                    }
                }
            }
            if (reachedDestination < 0) return null;

            var path = new List<int>();
            int node = reachedDestination;
            while (node != startNodeIndex)
            {
                int slot = previousSlot[node];
                path.Add(edgeIndexes[slot]);
                node = _outgoingAdjacency.FromNode[slot];
            }
            path.Reverse();
            return path;
        }
    }
}
