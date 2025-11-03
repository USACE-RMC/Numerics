/**
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* **/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

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
    /// <param name="edgeIndex">Index of the edge.</param>
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
        /// Index of the edge, often used as an index to the edge source (e.g., road segment).
        /// </summary>
        public int Index = edgeIndex;
    }

    /// <summary>
    /// Dijkstra dynamic programming implementation for shortest path optimization.
    /// </summary>
    public static class Dijkstra
    {
        private const int NEXT_NODE = 0;
        private const int EDGE_INDEX = 1;
        private const int COST = 2;

        /// <summary>
        /// May be a useful call in LifeSim -> GetPath().
        /// Follows the logic that is implemented in the Solve method.
        /// </summary>
        /// <param name="resultTable"></param>
        /// <param name="nodeIndex"></param>
        /// <returns></returns>
        public static bool PathExists(float[,] resultTable, int nodeIndex)
        {
            return !float.IsPositiveInfinity(resultTable[nodeIndex, COST]);
        }
        /// <summary>
        /// Solves the shortest path from every node in the network of edges to a given destination.
        /// </summary>
        /// <param name="edges">Edges, or segments, that make up the network.</param>
        /// <param name="destinationIndices">Indices of the destination nodes.</param>
        /// <param name="nodeCount">Optional number of nodes in the network. If not provided it will be calculated internally.</param>
        /// <param name="edgesFromNodes">Optional list of incoming edges from each node in the network. If not provided or mismatched with edges it will be calculated internally.</param>
        /// <returns>Lookup table of shortest paths from any given node.</returns>
        public static float[,] Solve(IList<Edge> edges, int[] destinationIndices, int nodeCount = -1, List<Edge>[] edgesFromNodes = null)
        {
            // Set optional parameters if required.
            int nNodes = (nodeCount == -1) ? (edges.Max(o => Math.Max(o.FromIndex,o.ToIndex)) + 1) : nodeCount;

            if (edgesFromNodes == null || edgesFromNodes.Length != nNodes)
            {
                edgesFromNodes = new List<Edge>[nNodes];
                //
                foreach(var edge in edges)
                {
                    edgesFromNodes[edge.ToIndex] ??= new List<Edge>();
                    edgesFromNodes[edge.ToIndex].Add(edge);
                }
            }


            float[,] resultTable = new float[nNodes, 3];
            //int[] nodeState = new int[nNodes]; //0 - Node hasn't been scanned yet, 1 - Node has been solved for, 2 - Node has been scanned into heap but not solved for.
            for(int i = 0; i < nNodes; i++)
            {
                resultTable[i, NEXT_NODE] = -1;
                resultTable[i, EDGE_INDEX] = -1;
                resultTable[i, COST] = float.PositiveInfinity;
            }

            for (int i = 0; i < destinationIndices.Length; i++)
            {
                int destinationIndex = destinationIndices[i];
                var partialResult = Solve(edges, destinationIndex, nNodes, edgesFromNodes);
                for(int j = 0; j < nNodes; j++)
                {
                    // Keep better path
                    if (partialResult[j, COST] < resultTable[j, COST])
                    {
                        resultTable[j, NEXT_NODE] = partialResult[j, NEXT_NODE];
                        resultTable[j, EDGE_INDEX] = partialResult[j, EDGE_INDEX];
                        resultTable[j,COST] = partialResult[j,COST];
                    }
                }
            }
            return resultTable;
           
        }

        /// <summary>
        /// Solves the shortest path from every node in the network of edges to a given destination.
        /// </summary>
        /// <param name="edges">Edges, or segments, that make up the network.</param>
        /// <param name="destinationIndex">Index of the destination node.</param>
        /// <param name="nodeCount">Optional number of nodes in the network. If not provided it will be calculated internally.</param>
        /// <param name="edgesToNodes">Optional list of incoming edges from each node in the network. If not provided or mismatched with edges it will be calculated internally.</param>
        /// <returns>Lookup table of shortest paths from any given node.</returns>
        public static float[,] Solve(IList<Edge> edges, int destinationIndex, int nodeCount = -1, List<Edge>[] edgesToNodes = null)
        {
            // Set optional parameters if required.
            int nNodes = (nodeCount == -1) ? (edges.Max(o => Math.Max(o.FromIndex, o.ToIndex)) + 1) : nodeCount;

            if (edgesToNodes == null || edgesToNodes.Length != nNodes)
            {
                edgesToNodes = new List<Edge>[nNodes];
                //
                foreach (var edge in edges)
                {
                    edgesToNodes[edge.ToIndex] ??= new List<Edge>();
                    edgesToNodes[edge.ToIndex].Add(edge);
                }
            }

            // Prepare results table with destination defined.
            float[,] resultTable = new float[nNodes, 3];
            int[] nodeState = new int[nNodes]; //0 - Node hasn't been scanned yet, 1 - Node has been solved for, 2 - Node has been scanned into heap but not solved for.
            float[] nodeWeightToDestination = new float[nNodes];

            //Initialize all nodes are unreachable
            for (int i = 0; i < nNodes; i++)
            {
                resultTable[i, NEXT_NODE] = -1;
                resultTable[i, EDGE_INDEX] = -1;
                resultTable[i, COST] = float.PositiveInfinity;
                nodeWeightToDestination[i] = float.PositiveInfinity;
            }

            BinaryHeap<Edge> heap = new BinaryHeap<Edge>(10000);

            resultTable[destinationIndex, NEXT_NODE] = destinationIndex; //Tail
            resultTable[destinationIndex, EDGE_INDEX] = -1; //edge index
            resultTable[destinationIndex, COST] = 0; //Cumulative Weight
            nodeWeightToDestination[destinationIndex] = 0;
            heap.Add(new BinaryHeap<Edge>.Node(0, destinationIndex, new Edge(destinationIndex, destinationIndex, 0, -1)));
            nodeState[destinationIndex] = 2;

           
            while (heap.Count > 0)
            {
                var node = heap.RemoveMin();
                int current = node.Index;
                float cost = node.Weight;

                if (nodeState[current] == 1)
                    continue;

                nodeState[current] = 1;

                //resultTable[current, NEXT_NODE] = node.Value.ToIndex;
                //resultTable[current, EDGE_INDEX] = node.Value.Index;
                //resultTable[current, COST] = cost;

                if (edgesToNodes[current] == null)
                    continue;

                foreach (var edge in edgesToNodes[current])
                {
                    int from = edge.FromIndex;
                    int to = edge.ToIndex;
                    float newCost = cost + edge.Weight;

                    if (newCost < nodeWeightToDestination[from])
                    {
                        nodeWeightToDestination[from] = newCost;
                        var newNode = new BinaryHeap<Edge>.Node(newCost, from, edge);

                        if (nodeState[from] != 2)
                        {
                            heap.Add(newNode);
                            nodeState[from] = 2;
                        }
                        else
                        {
                            heap.DecreaseKey(newNode);
                        }

                        resultTable[from, NEXT_NODE] = to;
                        resultTable[from, EDGE_INDEX] = edge.Index;
                        resultTable[from, COST] = newCost;
                      
                    }
                }
            }
                for (int i = 0; i < nNodes; i++)
                {
                    if (nodeState[i] == 0)
                        Console.WriteLine($"Node{i} is unreachable from destination {destinationIndex}");
                }
                return resultTable;
            }

        
    }
}
