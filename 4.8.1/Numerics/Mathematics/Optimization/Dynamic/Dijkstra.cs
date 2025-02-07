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
    public struct Edge
    {
        public int FromIndex;
        public int ToIndex;
        public float Weight;
        public int Index;

        public Edge(int fromNodeIndex, int toNodeIndex, float edgeWeight, int edgeIndex)
        {
            FromIndex = fromNodeIndex;
            ToIndex = toNodeIndex;
            Weight = edgeWeight;
            Index = edgeIndex;
        }
    }

    public class Networks
    {
        private readonly List<Edge>[] _outgoingEdges;
        private readonly List<Edge>[] _incomingEdges;
        private readonly int _nodeCount;
        private readonly int[] _destinationIndices;
        //private readonly RoadSegment[] _segments;
        private readonly Edge[] _edges;

        //public RoadSegment[] Segments { get => _segments; }
        public int[] DestinationIndices { get => _destinationIndices; }
        public List<Edge>[] IncomingEdges { get => _incomingEdges; }
        public List<Edge>[] OutgoingEdges { get => _outgoingEdges; }

        public static void ResultsToCSV(float[,] result, string outputFilePath)
        {
            System.IO.StreamWriter sw = new System.IO.StreamWriter(outputFilePath);
            //Headers
            sw.WriteLine($"From_Node, To_Node, Edge, Weight_To_Destination");
            //Data
            for (Int32 i = 0; i < result.GetLength(0); i++) { sw.WriteLine($"{i}, {result[i, 0]}, {result[i, 1]}, {result[i, 2]}"); }
            sw.Close(); sw.Dispose();
        }

        public static int[] PathEdges(float[,] result, int startNodeIndex)
        {
            List<int> edgeList = new List<int>();
            int nodeIndex = startNodeIndex;
            // generate the shortest path list from the results table
            while (result[nodeIndex, 2] != 0)
            {
                edgeList.Add((int)result[nodeIndex, 1]);
                nodeIndex = (int)result[nodeIndex, 0];
            }
            return edgeList.ToArray();
        }

        //public int[] PathRecords(float[,] result, int startNodeIndex)
        //{
        //    List<int> featureList = new List<int>();
        //    int nodeIndex = startNodeIndex;
        //    // generate the shortest path list from the results table
        //    while (result[nodeIndex, 2] != 0)
        //    {
        //        featureList.Add(_segments[(int)result[nodeIndex, 1]].RecordIndex);
        //        nodeIndex = (int)result[nodeIndex, 0];
        //    }
        //    return featureList.ToArray();
        //}

        //public int[] PathRecords(float[,] result, int startRecordIndex)
        //{
        //    List<int> featureList = new List<int>();
        //    for (int i = 0; i < result.GetLength(0); i++)
        //    {
        //        if (_segments[(int)result[i, 1]].RecordIndex == startRecordIndex)
        //        {
        //            int nodeIndex = i;
        //            // generate the shortest path list from the results table
        //            while (result[nodeIndex, 2] != 0)
        //            {
        //                featureList.Add(_segments[(int)result[nodeIndex, 1]].RecordIndex);
        //                nodeIndex = (int)result[nodeIndex, 0];
        //            }
        //            break;
        //        }
        //    }

        //    return featureList.ToArray();
        //}
    }


    public static class Dijkstra
    {
        /// <summary>
        /// Solves the shortest path from every node in the network of edges to a given destination.
        /// </summary>
        /// <param name="edges">Edges, or segments, that make up the network.</param>
        /// <param name="destinationIndices">Indices of the destination nodes.</param>
        /// <param name="nodeCount">Optional number of nodes in the network. If not provided it will be calculated internally.</param>
        /// <param name="edgesToNodes">Optional list of incoming edges from each node in the network. If not provided or mismatched with edges it will be calculated internally.</param>
        /// <returns>Lookup table of shortest paths from any given node.</returns>
        public static float[,] Solve(IList<Edge> edges, int[] destinationIndices, int nodeCount = -1, List<Edge>[] edgesToNodes = null)
        {
            // Set optional parameters if required.
            int nNodes = (nodeCount == -1) ? (edges.Max(o => o.FromIndex > o.ToIndex ? o.FromIndex : o.ToIndex) + 1) : nodeCount;

            if (edgesToNodes == null || edgesToNodes.Length != edges.Count)
            {
                edgesToNodes = new List<Edge>[nNodes];
                //
                for (int i = 0; i < edges.Count; i++)
                {
                    if (edgesToNodes[edges[i].ToIndex] == null) { edgesToNodes[edges[i].ToIndex] = new List<Edge>(); }
                    edgesToNodes[edges[i].ToIndex].Add(edges[i]);
                }
            }

            float[,] resultTable = new float[nNodes, 3];
            int[] nodeState = new int[nNodes]; //0 - Node hasn't been scanned yet, 1 - Node has been solved for, 2 - Node has been scanned into heap but not solved for.

            // Identify first valid destination index.
            int startIndex = -1;
            int destinationIndex = -1;
            for (int i = 0; i < destinationIndices.Length; i++)
            {
                if (edgesToNodes[destinationIndices[i]] != null && edgesToNodes[destinationIndices[i]].Count != 0)
                {
                    startIndex = i;
                    destinationIndex = destinationIndices[startIndex];
                    break;
                }
            }

            // Calculate shortest path for the first destination
            nodeState[destinationIndex] = 1;
            resultTable[destinationIndex, 0] = destinationIndex; //Tail
            resultTable[destinationIndex, 1] = -1; //edge index
            resultTable[destinationIndex, 2] = 0; //Cumulative Weight

            int previousValue = destinationIndex;
            if (edgesToNodes[previousValue] == null) { return resultTable; }

            //Add the nodes to the heap connected to the given destination node.
            int nodeIndex = 0;
            BinaryHeap<Edge> heap = new BinaryHeap<Edge>(10000);
            float[] nodeWeightToDestination = new float[nNodes];
            foreach (Edge edge in edgesToNodes[previousValue])
            {
                nodeIndex = edge.FromIndex;
                //Update the nodes current state.
                if (nodeState[nodeIndex] == 0)
                {
                    heap.Add(new BinaryHeap<Edge>.Node(edge.Weight, nodeIndex, edge));
                    nodeState[nodeIndex] = 2;
                    nodeWeightToDestination[nodeIndex] = edge.Weight;
                }
                else if (nodeState[nodeIndex] == 2)
                {
                    if (nodeWeightToDestination[nodeIndex] > edge.Weight)
                    {
                        nodeWeightToDestination[nodeIndex] = edge.Weight;
                        heap.Replace(new BinaryHeap<Edge>.Node(edge.Weight, nodeIndex, edge));
                    }
                }
            }

            //Solve for all nodes.
            BinaryHeap<Edge>.Node resultNode;
            float cumulativeWeight = 0;
            do
            {
                resultNode = heap.RemoveMin();
                previousValue = resultNode.Index;
                nodeState[previousValue] = 1;
                nodeWeightToDestination[previousValue] = resultNode.Weight;

                //Destination is the first result.
                resultTable[resultNode.Index, 0] = resultNode.Value.ToIndex;
                resultTable[resultNode.Index, 1] = resultNode.Value.Index;
                resultTable[resultNode.Index, 2] = resultNode.Weight;

                //Get all connections to previous node and add the weight to the heap
                if (edgesToNodes[previousValue] == null) { continue; }
                foreach (Edge edge in edgesToNodes[previousValue])
                {
                    nodeIndex = edge.FromIndex; //(previousValue == item.fromIndex) ? item.toIndex : item.fromIndex;
                                                //Update the nodes current state.
                    if (nodeState[nodeIndex] == 0)
                    {
                        cumulativeWeight = edge.Weight + resultNode.Weight;
                        heap.Add(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
                        nodeState[nodeIndex] = 2;
                        nodeWeightToDestination[nodeIndex] = cumulativeWeight;
                    }
                    else if (nodeState[nodeIndex] == 2)
                    {
                        cumulativeWeight = edge.Weight + resultNode.Weight;
                        if (nodeWeightToDestination[nodeIndex] > cumulativeWeight)
                        {
                            nodeWeightToDestination[nodeIndex] = cumulativeWeight;
                            heap.Replace(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
                        }
                    }
                }

            } while (heap.Count > 0);


            // Now calculate for all other destinations, time savings are reached by stopping the search vein once 
            // a distance becomes greater than a previously calculated distance.
            for (int i = startIndex + 1; i < destinationIndices.Length; i++)
            {
                heap = new BinaryHeap<Edge>(10000);
                nodeState = new int[nNodes];
                nodeWeightToDestination = new float[nNodes];

                destinationIndex = destinationIndices[i];
                nodeState[destinationIndex] = 1;
                resultTable[destinationIndex, 0] = destinationIndex; //Tail
                resultTable[destinationIndex, 1] = -1; //edge index
                resultTable[destinationIndex, 2] = 0; //Cumulative Weight

                if (edgesToNodes[previousValue] == null) { continue; }

                //Add the nodes to the heap connected to the given destination node.
                nodeIndex = 0;
                foreach (Edge edge in edgesToNodes[previousValue])
                {
                    nodeIndex = edge.FromIndex;
                    //Update the nodes current state.
                    if (nodeState[nodeIndex] == 0)
                    {
                        heap.Add(new BinaryHeap<Edge>.Node(edge.Weight, nodeIndex, edge));
                        nodeState[nodeIndex] = 2;
                        nodeWeightToDestination[nodeIndex] = edge.Weight;
                    }
                    else if (nodeState[nodeIndex] == 2)
                    {
                        if (nodeWeightToDestination[nodeIndex] > edge.Weight)
                        {
                            nodeWeightToDestination[nodeIndex] = edge.Weight;
                            heap.Replace(new BinaryHeap<Edge>.Node(edge.Weight, nodeIndex, edge));
                        }
                    }
                }

                //Solve for all nodes.
                do
                {
                    resultNode = heap.RemoveMin();
                    // If a faster route has already been found stop trying.
                    if (resultNode.Weight > resultTable[resultNode.Index, 2] && (resultTable[resultNode.Index, 0] > 0))
                    {
                        continue;
                    }

                    previousValue = resultNode.Index;
                    nodeState[previousValue] = 1;
                    nodeWeightToDestination[previousValue] = resultNode.Weight;

                    //Destination is the first result.
                    resultTable[previousValue, 0] = resultNode.Value.ToIndex;
                    resultTable[previousValue, 1] = resultNode.Value.Index;
                    resultTable[previousValue, 2] = resultNode.Weight;

                    //Get all connections to previous node and add the weight to the heap
                    if (edgesToNodes[previousValue] == null) { continue; }
                    foreach (Edge edge in edgesToNodes[previousValue])
                    {
                        nodeIndex = edge.FromIndex; //(previousValue == item.fromIndex) ? item.toIndex : item.fromIndex;
                                                    //Update the nodes current state.
                        if (nodeState[nodeIndex] == 0)
                        {
                            cumulativeWeight = edge.Weight + resultNode.Weight;
                            heap.Add(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
                            nodeState[nodeIndex] = 2;
                            nodeWeightToDestination[nodeIndex] = cumulativeWeight;
                        }
                        else if (nodeState[nodeIndex] == 2)
                        {
                            cumulativeWeight = edge.Weight + resultNode.Weight;
                            if (nodeWeightToDestination[nodeIndex] > cumulativeWeight)
                            {
                                nodeWeightToDestination[nodeIndex] = cumulativeWeight;
                                heap.Replace(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
                            }
                        }
                    }

                } while (heap.Count > 0);
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
            int nNodes = (nodeCount == -1) ? (edges.Max(o => o.FromIndex > o.ToIndex ? o.FromIndex : o.ToIndex) + 1) : nodeCount;

            if (edgesToNodes == null || edgesToNodes.Length != nodeCount)
            {
                edgesToNodes = new List<Edge>[nNodes];
                //
                for (int i = 0; i < edges.Count; i++)
                {
                    if (edgesToNodes[edges[i].ToIndex] == null) { edgesToNodes[edges[i].ToIndex] = new List<Edge>(); }
                    edgesToNodes[edges[i].ToIndex].Add(edges[i]);
                }
            }

            // Prepare results table with destination defined.
            float[,] resultTable = new float[nNodes, 3];
            int[] nodeState = new int[nNodes]; //0 - Node hasn't been scanned yet, 1 - Node has been solved for, 2 - Node has been scanned into heap but not solved for.

            nodeState[destinationIndex] = 1;
            resultTable[destinationIndex, 0] = destinationIndex; //Tail
            resultTable[destinationIndex, 1] = -1; //edge index
            resultTable[destinationIndex, 2] = 0; //Cumulative Weight

            int previousValue = destinationIndex;
            if (edgesToNodes[previousValue] == null) { return resultTable; }

            //Add the nodes to the heap connected to the given destination node.
            int nodeIndex = 0;
            BinaryHeap<Edge> heap = new BinaryHeap<Edge>(10000);
            float[] nodeWeightToDestination = new float[nNodes];
            foreach (Edge edge in edgesToNodes[previousValue])
            {
                nodeIndex = edge.FromIndex;
                //Update the nodes current state.
                if (nodeState[nodeIndex] == 0)
                {
                    heap.Add(new BinaryHeap<Edge>.Node(edge.Weight, nodeIndex, edge));
                    nodeState[nodeIndex] = 2;
                    nodeWeightToDestination[nodeIndex] = edge.Weight;
                }
                else if (nodeState[nodeIndex] == 2)
                {
                    if (nodeWeightToDestination[nodeIndex] > edge.Weight)
                    {
                        nodeWeightToDestination[nodeIndex] = edge.Weight;
                        heap.Replace(new BinaryHeap<Edge>.Node(edge.Weight, nodeIndex, edge));
                    }
                }
            }

            //Solve for all nodes.
            BinaryHeap<Edge>.Node resultNode;
            float cumulativeWeight = 0;
            do
            {
                resultNode = heap.RemoveMin();
                previousValue = resultNode.Index;
                nodeState[previousValue] = 1;
                nodeWeightToDestination[previousValue] = resultNode.Weight;

                //Destination is the first result.
                resultTable[previousValue, 0] = resultNode.Value.ToIndex;
                resultTable[previousValue, 1] = resultNode.Value.Index;
                resultTable[previousValue, 2] = resultNode.Weight;

                //Get all connections to previous node and add the weight to the heap
                if (edgesToNodes[previousValue] == null) { continue; }
                foreach (Edge edge in edgesToNodes[previousValue])
                {
                    nodeIndex = edge.FromIndex; //(previousValue == item.fromIndex) ? item.toIndex : item.fromIndex;
                                                //Update the nodes current state.
                    if (nodeState[nodeIndex] == 0)
                    {
                        cumulativeWeight = edge.Weight + resultNode.Weight;
                        heap.Add(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
                        nodeState[nodeIndex] = 2;
                        nodeWeightToDestination[nodeIndex] = cumulativeWeight;
                    }
                    else if (nodeState[nodeIndex] == 2)
                    {
                        cumulativeWeight = edge.Weight + resultNode.Weight;
                        if (nodeWeightToDestination[nodeIndex] > cumulativeWeight)
                        {
                            nodeWeightToDestination[nodeIndex] = cumulativeWeight;
                            heap.Replace(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
                        }
                    }
                }

            } while (heap.Count > 0);

            return resultTable;
        }


    }
}
