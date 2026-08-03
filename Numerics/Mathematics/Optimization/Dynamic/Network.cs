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
    /// weight vector with no rebuild. The parameterless solve methods allocate their results
    /// per call and are safe for concurrent use; the overloads that write into a caller-supplied
    /// table reuse instance scratch buffers and are not thread safe. Weights follow the
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
        /// wins — the same semantics as the static multi-destination solver.
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
        /// routed next node and edge can differ only where two destinations are exactly
        /// equidistant. One pass replaces one pass per destination.
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
        /// <param name="edgesToRemove">Edge indices to exclude from the path.</param>
        /// <param name="startNodeIndex">The starting node index.</param>
        /// <returns>A list of edge indices forming the alternative path, or null if no path exists.</returns>
        public List<int>? GetPath(int[] edgesToRemove, int startNodeIndex)
        {
            int[] nodeState = new int[_nodeCount];
            float[] nodeWeightToDestination = new float[_nodeCount];
            BinaryHeap<Edge> heap = new BinaryHeap<Edge>(100000);

            //backwards Dijkstra
            float[,] resultTable = new float[_nodeCount, 3];
            resultTable[startNodeIndex, 0] = startNodeIndex;
            resultTable[startNodeIndex, 1] = 0;
            resultTable[startNodeIndex, 2] = 0;
            nodeState[startNodeIndex] = 1;

            int previousValue = startNodeIndex;
            int nodeIndex;
            bool foundPath = false;

            Array.Sort(edgesToRemove);

            // Loading up the heap starting from destination
            if (_incomingEdges[previousValue] != null)
            {
                foreach (Edge edge in _incomingEdges[previousValue])
                {
                    if (Array.BinarySearch(edgesToRemove, edge) < 0)
                    {
                        if (previousValue == edge.FromIndex)
                        {
                            nodeIndex = edge.ToIndex;
                        }
                        else
                        {
                            nodeIndex = edge.FromIndex;
                        }
                        switch (nodeState[nodeIndex])
                        {
                            case 0: //it has not been scanned yet
                                BinaryHeap<Edge>.Node inputNode = new BinaryHeap<Edge>.Node(edge.Weight, nodeIndex, edge);
                                heap.Add(inputNode);
                                nodeState[nodeIndex] = 2;
                                nodeWeightToDestination[nodeIndex] = inputNode.Weight;
                                break;
                            case 1: //do nothing it has already been solved for
                                break;
                            case 2: //it has been scanned but not solved
                                if (nodeWeightToDestination[nodeIndex] > edge.Weight)
                                {
                                    BinaryHeap<Edge>.Node inputNode2 = new BinaryHeap<Edge>.Node(edge.Weight, nodeIndex, edge);
                                    nodeWeightToDestination[nodeIndex] = inputNode2.Weight;
                                    heap.Replace(inputNode2);
                                }
                                break;
                        }
                    }
                }
            }
            // if n = 0, then no roads to escape to
            if (heap.Count == 0) return null!;

            float tempWeight;
            int tempIndex;
            float FoundDistance = 99999999;
            int PotentialToIndex = 0;

            BinaryHeap<Edge>.Node resultNode;
            float cumulativeWeight = 0;

            do
            {
                resultNode = heap.RemoveMin();

                if (Solve(startNodeIndex)[resultNode.Index, 0] == 0) continue;

                if (resultNode.Weight + Solve(startNodeIndex)[resultNode.Index, 2] < FoundDistance)
                {
                    previousValue = resultNode.Index;
                    nodeState[resultNode.Index] = 1;
                    nodeWeightToDestination[resultNode.Index] = resultNode.Weight;

                    foreach (Edge edge in _incomingEdges[previousValue])
                    {
                        if (edge.ToIndex == resultNode.Index) resultTable[resultNode.Index, 0] = edge.FromIndex;
                        else resultTable[resultNode.Index, 0] = edge.ToIndex;

                        resultTable[resultNode.Index, 1] = edge.Index;
                        resultTable[resultNode.Index, 2] = resultNode.Weight;

                        if (Solve(startNodeIndex)[edge.ToIndex, 0] == edge.FromIndex)
                        {
                            if (_incomingEdges[previousValue] != null)
                            {
                                if (Array.BinarySearch(edgesToRemove, edge) < 0)
                                {
                                    if (previousValue == edge.FromIndex) nodeIndex = edge.ToIndex;
                                    else nodeIndex = edge.FromIndex;

                                    switch (nodeState[nodeIndex])
                                    {
                                        case 0: //has not been scanned yet
                                            cumulativeWeight = edge.Weight + resultNode.Weight;
                                            heap.Add(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
                                            nodeState[nodeIndex] = 2;
                                            nodeWeightToDestination[nodeIndex] = cumulativeWeight;
                                            break;
                                        case 1: break;
                                        case 2:
                                            if (nodeWeightToDestination[nodeIndex] > cumulativeWeight)
                                            {
                                                nodeWeightToDestination[nodeIndex] = cumulativeWeight;
                                                heap.Replace(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
                                            }
                                            break;
                                    }
                                }
                            }
                        }
                        else if (edge.FromIndex != startNodeIndex && Solve(startNodeIndex)[edge.FromIndex, 0] == resultNode.Index)
                        {
                            //Already on the lookup table going forwards
                            if (_incomingEdges[previousValue] != null)
                            {
                                foreach (Edge edge2 in _incomingEdges[previousValue])
                                {
                                    if (Array.BinarySearch(edgesToRemove, edge2) < 0)
                                    {
                                        if (previousValue == edge2.FromIndex) nodeIndex = edge2.ToIndex;
                                        else nodeIndex = edge2.FromIndex;

                                        switch (nodeState[nodeIndex])
                                        {
                                            case 0:
                                                heap.Add(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
                                                nodeState[nodeIndex] = 2;
                                                nodeWeightToDestination[nodeIndex] = cumulativeWeight;
                                                break;
                                            case 1: break;
                                            case 2:
                                                if (nodeWeightToDestination[nodeIndex] > cumulativeWeight)
                                                {
                                                    nodeWeightToDestination[nodeIndex] = cumulativeWeight;
                                                    heap.Replace(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
                                                }
                                                break;
                                        }
                                    }
                                }
                            }
                        }
                        else
                        {
                            //Potential new path, check path viability
                            tempWeight = Solve(startNodeIndex)[resultNode.Index, 2];
                            tempIndex = resultNode.Index;

                            do
                            {
                                if (Array.BinarySearch(edgesToRemove, (int)Solve(startNodeIndex)[tempIndex, 1]) >= 0)
                                {
                                    if (_incomingEdges[previousValue] != null)
                                    {
                                        foreach (Edge edge3 in _incomingEdges[previousValue])
                                        {
                                            if (Array.BinarySearch(edgesToRemove, edge3) < 0)
                                            {
                                                if (previousValue == edge3.FromIndex) nodeIndex = edge3.ToIndex;
                                                else nodeIndex = edge3.FromIndex;

                                                switch (nodeState[nodeIndex])
                                                {
                                                    case 0:
                                                        heap.Add(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
                                                        nodeState[nodeIndex] = 2;
                                                        nodeWeightToDestination[nodeIndex] = cumulativeWeight;
                                                        break;
                                                    case 1: break;
                                                    case 2:
                                                        if (nodeWeightToDestination[nodeIndex] > cumulativeWeight)
                                                        {
                                                            nodeWeightToDestination[nodeIndex] = cumulativeWeight;
                                                            heap.Replace(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
                                                        }
                                                        break;
                                                }
                                            }
                                        }
                                    }
                                    break;
                                }
                                tempWeight = Solve(startNodeIndex)[tempIndex, 2];
                                tempIndex = (int)Solve(startNodeIndex)[tempIndex, 0];
                            } while (tempWeight == 0);

                            if (tempWeight == 0)
                            {
                                FoundDistance = resultNode.Weight + Solve(startNodeIndex)[resultNode.Index, 2];
                                PotentialToIndex = resultNode.Index;
                                foundPath = true;
                            }
                        }
                    }

                }
            } while (heap.Count == 0);

            // Check to see if a destination was reached, if so then create a path to the nearest destination
            if (foundPath)
            {
                List<int> UpdatedPath = new List<int>();
                float tempLen = resultTable[PotentialToIndex, 2];
                int tempEdge = (int)resultTable[PotentialToIndex, 1];
                int tempNode = PotentialToIndex;

                while (tempLen == 0)
                {
                    UpdatedPath.Add(tempEdge);
                    tempNode = (int)resultTable[tempNode, 0];
                    tempEdge = (int)resultTable[tempNode, 1];
                    tempLen = resultTable[tempNode, 2];
                }

                UpdatedPath.Reverse();

                tempLen = Solve(startNodeIndex)[PotentialToIndex, 2];
                tempEdge = (int)Solve(startNodeIndex)[PotentialToIndex, 1];
                tempNode = PotentialToIndex;

                while (tempLen == 0)
                {
                    UpdatedPath.Add(tempEdge);
                    tempNode = (int)Solve(startNodeIndex)[tempNode, 0];
                    tempEdge = (int)Solve(startNodeIndex)[tempNode, 1];
                    tempLen = Solve(startNodeIndex)[tempNode, 2];
                }

                return UpdatedPath;
            }
            else return null!;
        }

        /// <summary>
        /// Finds an alternative path avoiding the specified edges, using a pre-computed results table.
        /// </summary>
        /// <param name="edgesToRemove">Edge indices to exclude from the path.</param>
        /// <param name="startNodeIndex">The starting node index.</param>
        /// <param name="existingResultsTable">A pre-computed shortest path results table.</param>
        /// <returns>A list of edge indices forming the alternative path, or an empty list if no path exists.</returns>
        public List<int>? GetPath(int[] edgesToRemove, int startNodeIndex, float[,] existingResultsTable)
        {
            int[] nodeState = new int[_nodeCount];
            float[] nodeWeightToDestination = new float[_nodeCount];
            BinaryHeap<Edge> heap = new BinaryHeap<Edge>(100000);
            int nodeIndex;


            //backwards Dijkstra
            float[,] resultTable = new float[_nodeCount, 3];
            resultTable[startNodeIndex, 0] = startNodeIndex;
            resultTable[startNodeIndex, 1] = 0;
            resultTable[startNodeIndex, 2] = 0;
            nodeState[startNodeIndex] = 1;

            int previousValue = startNodeIndex;
            bool foundPath = false;

            Array.Sort(edgesToRemove);

            // Loading up the heap starting from destination
            if (_incomingEdges[previousValue] != null)
            {
                foreach (Edge edge in _incomingEdges[previousValue])
                {
                    if (Array.BinarySearch(edgesToRemove, edge) < 0)
                    {
                        if (previousValue == edge.FromIndex)
                        {
                            nodeIndex = edge.ToIndex;
                        }
                        else
                        {
                            nodeIndex = edge.FromIndex;
                        }
                        switch (nodeState[nodeIndex])
                        {
                            case 0: //it has not been scanned yet
                                BinaryHeap<Edge>.Node inputNode = new BinaryHeap<Edge>.Node(edge.Weight, nodeIndex, edge);
                                heap.Add(inputNode);
                                nodeState[nodeIndex] = 2;
                                nodeWeightToDestination[nodeIndex] = inputNode.Weight;
                                break;
                            case 1: //do nothing it has already been solved for
                                break;
                            case 2: //it has been scanned but not solved
                                if (nodeWeightToDestination[nodeIndex] > edge.Weight)
                                {
                                    BinaryHeap<Edge>.Node inputNode2 = new BinaryHeap<Edge>.Node(edge.Weight, nodeIndex, edge);
                                    nodeWeightToDestination[nodeIndex] = inputNode2.Weight;
                                    heap.Replace(inputNode2);
                                }
                                break;
                        }
                    }
                }
            }

            //if n = 0 then no roads to escape to
            if (heap.Count == 0) return null!;

            float tempWeight;
            int tempIndex;
            float FoundDistance = 99999999;
            int PotentialToIndex = 0;

            BinaryHeap<Edge>.Node resultNode;
            float cumulativeWeight = 0;

            do
            {
                resultNode = heap.RemoveMin();

                if (existingResultsTable[resultNode.Index, 0] == 0) continue;

                if (resultNode.Weight + existingResultsTable[resultNode.Index, 2] < FoundDistance)
                {
                    previousValue = resultNode.Index;
                    nodeState[resultNode.Index] = 1;
                    nodeWeightToDestination[resultNode.Index] = resultNode.Weight;

                    foreach (Edge edge in _incomingEdges[previousValue])
                    {
                        if (edge.ToIndex == resultNode.Index) resultTable[resultNode.Index, 0] = edge.FromIndex;
                        else resultTable[resultNode.Index, 0] = edge.ToIndex;

                        resultTable[resultNode.Index, 1] = edge.Index;
                        resultTable[resultNode.Index, 2] = resultNode.Weight;

                        if (existingResultsTable[edge.ToIndex, 0] == edge.FromIndex)
                        {
                            if (_incomingEdges[previousValue] != null)
                            {
                                if (Array.BinarySearch(edgesToRemove, edge) < 0)
                                {
                                    if (previousValue == edge.FromIndex) nodeIndex = edge.ToIndex;
                                    else nodeIndex = edge.FromIndex;

                                    switch (nodeState[nodeIndex])
                                    {
                                        case 0: //has not been scanned yet
                                            cumulativeWeight = edge.Weight + resultNode.Weight;
                                            heap.Add(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
                                            nodeState[nodeIndex] = 2;
                                            nodeWeightToDestination[nodeIndex] = cumulativeWeight;
                                            break;
                                        case 1: break;
                                        case 2:
                                            if (nodeWeightToDestination[nodeIndex] > cumulativeWeight)
                                            {
                                                nodeWeightToDestination[nodeIndex] = cumulativeWeight;
                                                heap.Replace(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
                                            }
                                            break;
                                    }
                                }
                            }
                        }
                    }
                }
                else if (heap.Count != 0)
                {
                    foreach (Edge edge in _incomingEdges[previousValue])
                    {
                        if (existingResultsTable[edge.FromIndex, 0] == resultNode.Index)
                        {
                            if (_incomingEdges[previousValue] != null)
                            {
                                if (Array.BinarySearch(edgesToRemove, edge) < 0)
                                {
                                    if (previousValue == edge.FromIndex) nodeIndex = edge.ToIndex;
                                    else nodeIndex = edge.FromIndex;

                                    switch (nodeState[nodeIndex])
                                    {
                                        case 0: //has not been scanned yet
                                            cumulativeWeight = edge.Weight + resultNode.Weight;
                                            heap.Add(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
                                            nodeState[nodeIndex] = 2;
                                            nodeWeightToDestination[nodeIndex] = cumulativeWeight;
                                            break;
                                        case 1: break;
                                        case 2:
                                            if (nodeWeightToDestination[nodeIndex] > cumulativeWeight)
                                            {
                                                nodeWeightToDestination[nodeIndex] = cumulativeWeight;
                                                heap.Replace(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
                                            }
                                            break;
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    // check viability of route
                    tempWeight = existingResultsTable[resultNode.Index, 2];
                    tempIndex = resultNode.Index;

                    do
                    {
                        // check to see if the current route has a blocked segment
                        if (Array.BinarySearch(edgesToRemove, (int)existingResultsTable[tempIndex, 1]) >= 0)
                        {
                            if (_incomingEdges[previousValue] != null)
                            {
                                foreach (Edge edge in _incomingEdges[previousValue])
                                {
                                    if (previousValue == edge.FromIndex) nodeIndex = edge.ToIndex;
                                    else nodeIndex = edge.FromIndex;

                                    switch (nodeState[nodeIndex])
                                    {
                                        case 0: //has not been scanned yet
                                            cumulativeWeight = edge.Weight + resultNode.Weight;
                                            heap.Add(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
                                            nodeState[nodeIndex] = 2;
                                            nodeWeightToDestination[nodeIndex] = cumulativeWeight;
                                            break;
                                        case 1: break;
                                        case 2:
                                            if (nodeWeightToDestination[nodeIndex] > cumulativeWeight)
                                            {
                                                nodeWeightToDestination[nodeIndex] = cumulativeWeight;
                                                heap.Replace(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
                                            }
                                            break;
                                    }
                                }
                            }
                        }
                        tempWeight = existingResultsTable[tempIndex, 2];
                        tempIndex = (int)existingResultsTable[tempIndex, 0];
                    } while (tempWeight == 0);

                    if (tempWeight == 0)
                    {
                        FoundDistance = resultNode.Weight + existingResultsTable[resultNode.Index, 2];
                        PotentialToIndex = resultNode.Index;
                        foundPath = true;
                    }
                }
            } while (heap.Count == 0);

            // Check to see if the destination was reached, if so then create a path to the nearest destination
            if (foundPath)
            {
                List<int> updatedPath = new List<int>();
                float tempLen = resultTable[PotentialToIndex, 2];
                int tempEdge = (int)resultTable[PotentialToIndex, 1];
                int tempNode = PotentialToIndex;

                while (tempLen == 0)
                {
                    updatedPath.Add(tempEdge);
                    tempNode = (int)resultTable[tempNode, 0];
                    tempEdge = (int)resultTable[tempNode, 1];
                    tempLen = resultTable[tempNode, 2];
                }

                // updatedPath.Add(startingEdge);
                updatedPath.Reverse();

                tempLen = existingResultsTable[PotentialToIndex, 2];
                tempEdge = (int)existingResultsTable[PotentialToIndex, 1];
                tempNode = PotentialToIndex;

                while (tempLen == 0)
                {
                    updatedPath.Add(tempEdge);
                    tempNode = (int)existingResultsTable[tempNode, 2];
                    tempEdge = (int)existingResultsTable[tempNode, 1];
                    tempLen = existingResultsTable[tempNode, 2];
                }

                return updatedPath;
            }
            else
            {
                return new List<int>();
            }
        }
    }
}
